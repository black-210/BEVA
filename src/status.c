#include "beva/status.h"
const char *beva_status_string(beva_status s) {
    switch(s) {
        case BEVA_OK: return "ok";
        case BEVA_ERR_NULL: return "null argument";
        case BEVA_ERR_ALLOC: return "allocation failure";
        case BEVA_ERR_DIMENSION: return "dimension mismatch";
        case BEVA_ERR_SINGULAR: return "singular matrix";
        case BEVA_ERR_INVALID: return "invalid argument";
        case BEVA_ERR_OVERFLOW: return "integer overflow";
        case BEVA_ERR_CONVERGENCE: return "failed to converge";
        case BEVA_ERR_DOMAIN: return "domain error";
        case BEVA_ERR_NOT_IMPLEMENTED: return "not implemented";
        default: return "unknown error";
    }
}
