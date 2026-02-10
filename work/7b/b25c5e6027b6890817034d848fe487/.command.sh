#!/bin/bash -ue
set -euo pipefail
echo "✗ No index found for subset-recalc-pop12_sub1_chr21.vcf.gz — creating..." > "subset-recalc-pop12_sub1_chr21.log"
tabix -p vcf "subset-recalc-pop12_sub1_chr21.vcf.gz"
nvariant="$(bcftools index -n subset-recalc-pop12_sub1_chr21.vcf.gz)"
nsample="$(bcftools query -l subset-recalc-pop12_sub1_chr21.vcf.gz | wc -l)"
ngenotypes="$(( nsample * nvariant ))"
{
  echo " === EGA BCFTools Pipeline ==="
  echo "Developed by: Mireia Marin Ginestar (mireia.marin@crg.eu)"
  echo "version 3.0.0"
  echo ""
  echo "✓ Index created for subset-recalc-pop12_sub1_chr21.vcf.gz"
  echo ""
  echo " === ORIGINAL STATISTICS === "
  echo "Variant Number: $nvariant"
  echo "Sample Number: $nsample"
  echo "Genotype Number: $ngenotypes"
} >> "subset-recalc-pop12_sub1_chr21.log"
