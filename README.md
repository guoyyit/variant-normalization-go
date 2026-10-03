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
| `cmd/varnorm` | `import` / `resume` / `query` 命令行 |

## 规范化规则

规范化实现 Tan 等人 2015（vt）/ bcftools `norm` 的修剪与左对齐算法：

1. **校验**：`CHROM` 非空、`POS ≥ 1`、REF/ALT 仅含 `ACGTN`（大写归一）、
   ALT 与 REF 不同、基因型索引在范围内。
2. **坐标越界检查**：当提供参考基因组时，`POS + len(REF) − 1` 超过染色体长度
   即拒绝。
3. **符号等位**：`<DEL>`、breakend（含 `[` `]`）、`*` 只做坐标校验并原样保留。
4. **修剪**：先去共同后缀（保留至少一个锚定碱基），再去共同前缀（锚定碱基
   紧邻变异左侧）。
5. **左对齐**：对 indel 逐碱基把等位窗口向左移动，直到无法继续（撞到染色体起点、
   两条等位等长，或左端碱基不再相同）。每步向左端补一个参考碱基（保持单倍型
   不变），再修剪共同后缀。结果满足幂等：对规范结果再次规范化不变。

### 规范化等价性定义

规范化必须保持语义等价：对任意合法输入，把输入记录与输出记录分别应用到同一
参考序列上重建单倍型序列，两条序列必须完全一致。规范化只允许改变表示形式，
不允许改变变异本身。

同一变异的所有等价文本形式必须收敛到唯一规范表示，得到完全相同的
`(CHROM, POS, REF, ALT)`。对规范结果再次规范化，结果保持不变（幂等）。

等价性验收独立于规范化实现：自动测试把输入与输出分别应用到参考序列上重建
单倍型并比较，覆盖程序化生成的枚举空间（有限字母表上所有长度受限的等位基因对
与位置组合）。

同一变异无论以何种文本形式上报，规范化后都得到相同的
`(CHROM, POS, REF, ALT)`。错误通过 `normalize.Error` 区分：

- `errors.Is(err, normalize.ErrInvalidFormat)` — 格式非法；
- `errors.Is(err, normalize.ErrOutOfBounds)` — 坐标越界；

`Error.Reason` 给出人类可读原因（非法字段、越界长度等）。

## 存量数据重建

修复前生成的存量库中的记录按旧规范化逻辑写入。修复后继续导入或重导同一批
数据时，同一变异不得以不同规范表示出现两条记录。需要提供存量数据重建方式，
使存量库获得修复后的规范表示。

`varnorm rebuild` 重新导入原始 VCF 文件到全新存储，再原子替换存量库：

1. 在临时目录（`<db>.rebuild`）创建新存储，保持存量库可正常查询。
2. 流式重新导入原始 VCF 文件（内存占用受 memtable 限制，与总量无关）。
3. 进度按 `(source, file)` 持久化到检查点，崩溃后从最近检查点继续。
4. 导入完成后关闭新存储，原子替换存量库（rename）。

重建期间查询持续可用：查询存量库看到旧完整状态，替换后看到新完整状态，
不会看到部分重建的中间状态。重建是显式操作，打开存量库时不会静默改写数据。

重建完成后存量库中每条记录与用修复后逻辑全新导入同一批原始数据得到的记录
完全一致。

## 验证方法

- **语义等价性枚举**：`go test ./pkg/normalize/ -run TestEquivalenceEnumeration`
  枚举有限字母表上所有长度受限的等位基因对与位置组合，把输入与输出分别应用
  到参考序列重建单倍型并比较。
- **同变异多形式收敛**：`go test ./pkg/normalize/ -run TestCanonical`
  验证同一变异的不同文本形式收敛到唯一规范表示。
- **跨来源合并与区间查询端到端**：`go test ./pkg/importer/ -run TestE2E`
  验证跨来源合并、区间查询正确性、排序与翻页稳定性。
- **存量重建与中断恢复**：`go test ./pkg/importer/ -run TestRebuild`
  验证重建结果与全新导入一致、崩溃恢复不重复不遗漏、重建与查询并发。
- **并发安全**：`go test -race ./...` 验证写入与查询并发时无数据竞争。

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

# 存量数据重建（重新导入原始 VCF，原子替换存量库）
varnorm rebuild -source lab-a -reference ref.fa -config config.json \
  -db varnorm.db -checkpoint cp.json input.vcf [input2.vcf ...]
```

非法行不会中断导入：它们计入报告的 `rejected` 并带行号与原因输出到 stderr。

## 测试

`go test ./...` 覆盖：规范化（SNP/MNP/插入/删除/左对齐/幂等/越界/非法）、
语义等价性枚举、同变异多形式收敛、多等位拆分与基因型重映射、跨来源去重合并、
区间查询正确性、排序与翻页稳定性、乱序到达、QUAL/FILTER 策略、批量检查点与
崩溃恢复、流式导入内存有界、存量重建与中断恢复，以及写入与区间查询并发（`-race`）。
