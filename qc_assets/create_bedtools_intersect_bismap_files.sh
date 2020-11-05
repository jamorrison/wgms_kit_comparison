# Find which positions in the bismap multi-read mapabbility file overlap with
# the different region files

mkdir -p bismap

bedtools intersect -a k100.bismap.bedgraph.gz -b hg38/cgi.bed.gz |
gzip -c > bismap/cgi_bismap.bed.gz
bedtools intersect -a k100.bismap.bedgraph.gz -b hg38/cpg.bed.gz |
gzip -c > bismap/cpg_bismap.bed.gz
bedtools intersect -a k100.bismap.bedgraph.gz -b hg38/exon.bed.gz |
gzip -c > bismap/exon_bismap.bed.gz
bedtools intersect -a k100.bismap.bedgraph.gz -b hg38/genic_regions.bed.gz |
gzip -c > bismap/genic_regions_bismap.bed.gz
bedtools intersect -a k100.bismap.bedgraph.gz -b hg38/intergenic_regions.bed.gz |
gzip -c > bismap/intergenic_regions_bismap.bed.gz
bedtools intersect -a k100.bismap.bedgraph.gz -b hg38/rmsk.bed.gz |
gzip -c > bismap/rmsk_bismap.bed.gz
