# varnorm — VCF 变异规范化与多来源合并

`varnorm` 接收多个来源上报的 VCF 文本记录，把它们规范化为唯一表示、拆分多等位
记录、跨来源去重合并、重算注释与统计字段，并支持按染色体区间稳定分页查询，以及
可断点恢复的流式批量导入。

## 架构

| 包 | 职责 |
| --- | --- |
| `pkg/vcf` | 流式 VCF 解析（逐行，固定缓冲），输出 `variant.Record` |
| `pkg/normalize` | 合法性校验、碱基大写、等位修剪、indel 左对齐；区分两类拒绝原因 |
| `pkg/split` | 多等位记录拆分为单等位记录，基因型与 FORMAT 字段重映射 |
| `pkg/genome` | FASTA 参考序列加载（左对齐查询） |
| `pkg/merge` | 顺序无关、确定性的记录合并与来源追溯 |
| `pkg/stats` | AC/AN/AF、DP、QUAL、FILTER 重算与阈值过滤 |
| `pkg/config` | 来源优先级、QUAL 聚合、过滤阈值等 JSON 配置 |
| `pkg/store` | 内存表 + WAL + 磁盘段 + 键索引 + 后台压实的持久化存储 |
| `pkg/importer` | 解析→规范化→拆分→原子合并入库→批量检查点的导入管线 |
| `pkg/rebuild` | 用原始 VCF 重建影子世代并原子切换的存量重建 |
| `cmd/varnorm` | `import` / `resume` / `query` 命令行 |

## 规范化规则

规范化实现 Tan 等人 2015（vt）/ bcftools `norm` 的修剪与左对齐算法：

1. **校验**：`CHROM` 非空、`POS ≥ 1`、REF/ALT 仅含 `ACGTN`（大写归一）、
   ALT 与 REF 不同、基因型索引在范围内。
2. **坐标越界检查**：当提供参考基因组时，`POS + len(REF) − 1` 超过染色体长度
   即拒绝。
3. **符号等位**：`<DEL>`、breakend（含 `[` `]`）、`*` 只做坐标校验并原样保留。
4. **修剪**：先去共同后缀（POS 不变），再去共同前缀（每去掉一个共同前缀碱基，
   POS 向右推进 1），始终保留至少一个锚定碱基。
5. **左对齐**：对 indel 逐碱基把等位窗口向左移动。窗口左移一等价于
   “在左端加入共同锚定碱基、在右端各丢弃一个碱基”，因此只有当两条等位**右端被
   丢弃的碱基相同**（移位前 `REF` 与 `ALT` 末位碱基相等）时才能左移；每移一步都
   重新修剪共同前后缀。直到撞到染色体起点、两条等位等长、或右端碱基不一致为止。
   结果满足幂等：对规范结果再次规范化不变。

同一变异无论以何种文本形式上报，规范化后都得到相同的
`(CHROM, POS, REF, ALT)`。错误通过 `normalize.Error` 区分：

- `errors.Is(err, normalize.ErrInvalidFormat)` — 格式非法；
- `errors.Is(err, normalize.ErrOutOfBounds)` — 坐标越界；

`Error.Reason` 给出人类可读原因（非法字段、越界长度等）。

## 规范化的语义等价性定义

规范化**只允许改变表示形式，不允许改变变异本身**。对任意合法记录
`(POS, REF, ALT)`，把输入与规范结果分别应用到同一参考序列 `R` 上重建单倍型：

```
H = R[:POS-1] + ALT + R[POS-1+len(REF):]
```

两条单倍型必须逐碱基一致。等价的不同上报形式（indel 在重复序列中可有多个合法
对齐、共同锚定碱基可多可少）必须收敛到唯一的 `(CHROM, POS, REF, ALT)`，且该结果
再次规范化保持不变（幂等）。

历史上左对齐无条件移动窗口、且前缀修剪不推进 POS，会把变异改写为另一个变异
（如 `ACGT` 上 pos2 的 `C>CA` 插入 A 被错误改成 pos1 的 `A>AC` 插入 C）。当前
实现严格按上面的等价条件左移，并用独立重建做验收。

### 等价性验收（独立于实现）

`pkg/normalize/equiv_test.go` 不使用被测实现的任何中间结果作为期望值：

- 用独立的 `reconstructHaplotype` 分别对输入与输出重建单倍型并逐碱基比较；
- 程序化枚举有限字母表（`{A,C}`）上长度 9 的全部参考序列、长度 1–3 的全部
  REF/ALT 等位与全部合法位置（约 16 万组合）；
- 对每个变体，独立枚举参考单倍型与重建单倍型之间所有合法的公共前缀/后缀切分，
  得到该等价类的全部“已修剪”拼写，并独立选取最小 `(POS, REF, ALT)` 作为规范键；
- 断言：输入/输出单倍型一致、规范化幂等、等价类全部拼写收敛且等于独立最小键。

运行：`go test ./pkg/normalize/ -run HaplotypeEquivalence -v`。

## 多等位拆分

`split.Split` 把一条含 N 个 ALT 的记录拆成 N 条单等位记录：

- 第 `i` 条拆分的基因型把原 ALT 索引 `i+1` 重映射为 `1`，REF 保持 `0`，其它
  等位重映射为缺失 `-1`；相位（phased）保留。
- FORMAT 字段按 VCF 基数约定尽力重映射：`AD` 这类含 REF 的 R 型字段保留
  `REF,选中ALT`；其它逗号分隔的 A 型字段取选中 ALT 的值；标量字段原样保留。
- 拆分前的 `AC/AN/AF` 被删除，由合并后的重算步骤重新生成，避免遗留矛盾值。
- 拆分前后等位基因总数不变（信息保持）。

## 合并策略

合并以规范化后的 `(CHROM, POS, REF)` 为分组键，同组记录调用 `merge.Merge`，
结果**与到达顺序无关且可解释**：

- **ALT 槽位**：取两侧 ALT 并集并按字典序排序赋槽；相同 ALT 合并到同一槽。
- **基因型**：样本按名归并，等位索引按新槽位重编号；同一相位槽位的冲突调用
  值由**来源优先级**高的记录胜出，优先级相同按稳定顺序裁决。
- **QUAL**：按 `qual_agg` 聚合——`max`（默认）、`min`、`weighted`
  （用两侧 `INFO/DP` 加权；缺权重时退化为 max）。
- **FILTER**：见下方 `filter_mode`。
- **INFO**：等位无关的标量键在两侧一致时保留，仅一侧提供时保留；冲突值不
  静默选择而是删除；`AC/AN/AF` 一律由 `stats.Recompute` 重新计算。
- **ID**：两侧 `ID` 去重并字典序排序。
- **来源追溯**：保留全部不同的 `(source, file)` 来源；同一条来源对同一键
  的重复上报在折叠函数中短路，保证重放/断点续传幂等。

合并后 `stats.Recompute` 基于实际样本基因型重算：

- `AN` = 非缺失染色体数；`AC[i]` = ALT 槽 `i` 的计数；`AF[i] = AC[i]/AN`；
  `DP` = 样本 `FORMAT/DP` 之和。
- `min_qual`：`QUAL < min_qual` 时加 `MinQual` 过滤码；
  `min_af`：最大 AF 低于阈值时加 `MinAF` 过滤码；无任何过滤码时为 `PASS`。
  历史的、与合并后内容矛盾的 QUAL/FILTER/AC/AN/AF 不会被保留。

## 配置项

见 `config.example.json`：

| 键 | 含义 |
| --- | --- |
| `id` | 配置版本号，写入每条合并记录的 `ConfigID` 用于追溯 |
| `source_priority` | 来源名 → 优先级（整数，越大越高） |
| `default_priority` | 未知来源的默认优先级 |
| `min_qual` | QUAL 阈值，低于则 `MinQual` |
| `min_af` | AF 阈值（`[0,1]`），最高 AF 低于则 `MinAF` |
| `qual_agg` | `max`（默认）/ `min` / `weighted` |
| `filter_mode` | `all`（默认，全部来源 PASS 才 PASS）/ `any`（任一来源 PASS 即 PASS） |

配置只在导入器构造时读取并固化到折叠函数；**配置变更只影响之后的导入**，
已入库记录保持自洽（其 `ConfigID` 记录生成时所用配置）。

## 查询

`store.Query` 按染色体与 1-based 闭区间 `[start,end]` 查询，结果按
`(chrom, pos, ref)` 稳定排序。分页使用不透明 base64 token（本质是上一页末键），
翻页稳定。写路径在写锁内原子完成「读已有→合并→写 WAL/内存表」，因此并发查询
不会观察到合并的中间状态；读路径持有 RLock，与写入互斥但与其它读并发。

## 存储与流式导入

- 内存使用受 `WithMemLimit`（内存表键数）限制，与输入总量无关：内存表满后
  原子刷成不可变磁盘段，后台在段数累积时归并压实。
- 键索引只保存 `key → 段指针`（不缓存记录体），内存只随**去重后键数**增长；
  记录体按需读盘。
- WAL 采用组提交：导入器按批（`WithSyncBatch`，默认 1000 条）先 `Sync`
  WAL 再推进检查点，避免逐条 fsync；正常关闭时强制落盘。
- 检查点按来源保存 `(字节偏移, 行号)`，原子临时文件 + rename + fsync 持久化，
  崩溃只可能留下旧或新完整清单。`resume` 从最近检查点继续；由于折叠按
  `(source,file,line)` 去重，一个批次的重放不重复、不遗漏。

## 存量数据重建

修复规范化语义后，旧版本写入的库可能保存了被改写的变异。由于被改写的记录已经
丢失原始含义，无法只靠库内数据还原，因此重建被定义为**用原始上报 VCF 重新导入
一个全新的影子世代**：这正是“与用修复后逻辑全新导入同一批原始数据完全一致”的
构造方式。

- 数据库根目录下的 `current.json` 是世代指针；没有该文件的旧库按根目录直接打开，
  **打开、查询、普通导入都不会静默改写旧数据**。
- `varnorm rebuild` 写入新的 `gen-<序号>/` 影子世代，活动世代在整个过程中持续
  提供查询；重建完成时先 fsync 影子世代，再原子 rename `current.json`，最后在锁
  保护与原子指针交换下切换活动库。查询只能观察到重建前或重建后的**完整**状态。
- 重建是流式的：内存受影子库内存表上限约束，不随记录总量线性增长。
- 重建进度持久化到独立检查点及其 `.rebuild` 影子世代状态文件；中途失败后再次
  `rebuild` 会复用同一影子世代并从最近检查点继续，折叠按来源去重，重放不重复、
  不遗漏。已提交的旧世代文件保留为备份，只是不再被引用。
- 重建是**显式操作**：`import`/`resume`/`query` 的调用方式与配置格式不变，未
  重建的旧库仍可正常打开与查询。

验证重建结果：对同一批原始文件，分别执行一次“全新导入临时库”和“对存量库
rebuild”，两个库逐条记录比较应完全一致（`pkg/rebuild` 的端到端测试即按此断言）。

## 命令行

```sh
# 导入（需要左对齐时提供 FASTA 参考）
varnorm import -source lab-a -reference ref.fa -config config.json \
  -db varnorm.db -checkpoint cp.json input.vcf

# 中断后从检查点继续
varnorm resume -source lab-a -reference ref.fa -config config.json \
  -db varnorm.db -checkpoint cp.json input.vcf

# 区间查询，稳定分页
varnorm query -db varnorm.db -chrom chr1 -start 100 -end 200 -page-size 50
varnorm query -db varnorm.db -chrom chr1 -start 100 -end 200 \
  -page-token "<token>" -json

# 显式存量重建（修复规范化后，用原始 VCF 重建新世代）
varnorm rebuild -source lab-a -reference ref.fa -config config.json \
  -db varnorm.db -checkpoint rebuild-cp.json input1.vcf [input2.vcf ...]
# 中断后重复执行同一条命令即可从检查点继续，直到发布新世代
```

非法行不会中断导入：它们计入报告的 `rejected` 并带行号与原因输出到 stderr。

## 测试

`go test ./...` 覆盖：规范化（SNP/MNP/插入/删除/左对齐/幂等/越界/非法）、
多等位拆分与基因型重映射、跨来源去重合并、乱序到达、QUAL/FILTER 策略、
批量检查点与崩溃恢复、流式导入内存有界，以及写入与区间查询并发（`-race`）；
此外还包括：约 16 万组合的**单倍型等价性枚举/收敛/幂等**（期望值独立重建推导）、
同变异多形式合并与区间命中间隔、存量重建与全新导入逐记录一致、重建中断续传、
以及重建与查询并发时只能观察到完整世代。
