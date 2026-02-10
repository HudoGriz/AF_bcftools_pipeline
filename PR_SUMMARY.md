# Fix Log File Staging Issues

## Problem
The pipeline was experiencing Nextflow staging errors when trying to pass the same log file through multiple processes. The error occurred because Nextflow cannot stage the same file multiple times as both input and output in different processes.

## Solution
Refactored log file handling throughout the pipeline to eliminate staging conflicts:

### Key Changes

1. **Log File Publishing Strategy**
   - Added `publishDir` directive with `mode: 'copy'` to all processes that handle logs
   - Changed output patterns from `path(log)` to `path("*.log")` for glob matching
   - This ensures logs are properly published to the output directory at each step

2. **Log File Propagation**
   - Each process now creates its own copy of the log file from the previous stage
   - Added `cp ${log} ${base_name}.log` at the start of each script
   - Base name extraction using regex to remove process-specific suffixes

3. **Base Name Handling**
   - Added `def base_name` variable in each process to extract consistent log names
   - Uses regex patterns like `.replaceAll(/_split-multiallelic$/, '')` to strip suffixes
   - Ensures log files maintain their original base name throughout the pipeline

### Modified Processes

All processes were updated with consistent log handling:
- `INDEX_VCF` - Creates initial log file
- `SPLIT_MULTIALLELIC` - Copies and appends to log
- `MASK_GT` - Copies and appends to log  
- `VARIANT_QC` - Copies and appends to log
- `SAMPLE_QC` - Copies and appends to log
- `FIX_PLOIDY` - Copies and appends to log
- `ADD_AF` - Copies and appends to log (final output with mode 'move' for VCF)

### Benefits

✅ Eliminates Nextflow staging conflicts  
✅ Each process operates on its own log file copy  
✅ Logs are properly published at each step  
✅ Base names are preserved throughout the pipeline  
✅ All log entries accumulate in a single file per sample  
✅ Clean separation between intermediate files and final outputs

## Testing

The changes have been linted and pass all Nextflow syntax checks:
```bash
nextflow lint main.nf
# ✅ 1 file had no errors
```

## Files Modified

- `main.nf` - Updated all 7 processes with new log handling logic
- `nextflow.config` - Added Docker support and container configuration for testing

## Next Steps

After merging, the pipeline should run without log staging errors. Future improvements could include:
- Implementing nf-core modules for standardized bcftools operations
- Adding comprehensive test suite with CI/CD
- Container optimization and multi-platform support
