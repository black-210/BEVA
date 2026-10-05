#ifndef BEVA_STATUS_H
#define BEVA_STATUS_H
#ifdef __cplusplus
extern "C" {
#endif
typedef enum {
    BEVA_OK = 0,
    BEVA_ERR_NULL = 1,
    BEVA_ERR_ALLOC = 2,
    BEVA_ERR_DIMENSION = 3,
    BEVA_ERR_SINGULAR = 4,
    BEVA_ERR_INVALID = 5,
    BEVA_ERR_OVERFLOW = 6,
    BEVA_ERR_CONVERGENCE = 7,
    BEVA_ERR_DOMAIN = 8,
    BEVA_ERR_NOT_IMPLEMENTED = 9
} beva_status;
const char *beva_status_string(beva_status s);
#ifdef __cplusplus
}
#endif
#endif
