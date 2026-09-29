package store

import (
	"encoding/base64"
	"errors"
)

// encodeToken wraps an opaque sort key as a base64 page token.
func encodeToken(key string) string {
	return base64.RawURLEncoding.EncodeToString([]byte(key))
}

func decodeToken(token string) (string, error) {
	raw, err := base64.RawURLEncoding.DecodeString(token)
	if err != nil {
		return "", errors.New("invalid page token")
	}
	return string(raw), nil
}
