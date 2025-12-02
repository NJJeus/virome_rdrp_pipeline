

curl -o NeoRdRp.2.1.hmm.xz https://zenodo.org/records/10851672/files/NeoRdRp.2.1.hmm.xz\?download\=1
curl -o Annotation_of_Seed_RdRp_datasets.tsv.xz https://zenodo.org/records/10851672/files/Annotation_of_Seed_RdRp_datasets.tsv.xz\?download\=1

xz -d NeoRdRp.2.1.hmm.xz
xz -d Annotation_of_Seed_RdRp_datasets.tsv.xz

mkdir HMMs

csplit -z -f HMMs/RdRp NeoRdRp.2.1.hmm /^HMMER3/ '{*}'


