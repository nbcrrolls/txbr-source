include ../ConfigInfo.txt
include ../MakeConfig

INCLUDE_FLAGS := $(IMOD_INCLUDE_FLAGS) $(TXBR_INCLUDE_FLAGS)

all: deviceQuery

deviceQuery: deviceQuery.cu
	$(NVCC) deviceQuery.cu $(CUDA_INCLUDE_FLAGS) -o deviceQuery $(CUDA_LIB_FLAGS)

clean:
	find . -type f -name "*.o"  -exec rm "{}" \;
	find . -type f -name "*.so" -exec rm "{}" \;
	find . -type f -name "*.linkinfo" -exec rm "{}" \;
	find . -type f -name "*~"  -exec rm "{}" \;
	find . -type f -name "core*"  -exec rm "{}" \;
	find . -type f -name "deviceQuery" -exec rm "{}" \;

