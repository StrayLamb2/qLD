#define LOCAL_0 16
#define LOCAL_1 2

#define GPU_BLOCK_NC 8192
#define GPU_BLOCK_KC 2048
#define GPU_BLOCK_MC 8192

#define GPU_BLOCK_NR 128
#define GPU_BLOCK_MR 128

// #define BLOCK_SIZE_X 128
// #define BLOCK_SIZE_Y 128

#define BLOCK_SIZE_X 128
#define BLOCK_SIZE_Y 8

// NC/MC 32, n/m 32 and 64,64 work fine (k = 1)
// NC/MC 32, n/m 64 does not
