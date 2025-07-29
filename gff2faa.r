#!/usr/bin/env Rscript

input<-commandArgs(trailingOnly=TRUE)
f.ref<-input[1]
f.gff<-input[2]
f.out.gene<-input[3]
f.out.prot<-input[4]

if(is.na(f.ref))
{
	cat('
# Get protein-coding sequences of genes using reference genome and gff file.
#
# Usage: gff2faa.r $fasta $gff $out_DNA $out_protein
#	fasta: reference genome file in fasta format.
#	gff: gene annotation file in gff format.
#	out_DNA: output file for DNA sequences of the protein-coding genes.
#	out_protein: output file for protein sequences of the protein-coding genes (stop codon removed).
#
# Requirements for the gff file:
#	CDS tag for Column 3;
#	Parent tag for Column 9.
#
# Dependency: Biostrings.

')
	quit("no")
}

suppressPackageStartupMessages(library(Biostrings))

# Get protein-coding sequences of genes using reference genome and gff file.
# Requirements:
#	CDS tag for Column 3;
#	Parent tag for Column 9.
# gff: a data.frame of gff format.
# ref: a DNAStringSet of reference genome.
# Returns a DNAStringSet of protein-coding sequences of genes.
gff2fa<-function(gff,ref)
{
	names(gff)<-sprintf("C%d",1:ncol(gff))
	index<-gff$C1%in%names(ref)
	gff<-gff[index,]
	index<-gff$C3=="CDS"
	gff<-gff[index,]
	index<-with(gff,order(C1,C4,C5))
	gff<-gff[index,]
	
	cds<-with(gff,subseq(ref[C1],C4,C5))
	srd<-gff$C7
	cdn<-gff$C8
	parent<-sub("^.*Parent=([^;]*).*$","\\1",gff$C9)
	
	cds.list<-split(cds,parent)
	srd.list<-split(srd,parent)
	cdn.list<-split(cdn,parent)
	
	gene<-sapply(cds.list,paste,collapse="")
	gene<-DNAStringSet(gene)
	strd<-sapply(srd.list,function(x)x[1])
	
	index<-strd=="-"
	gene[index]<-reverseComplement(gene[index])
	cdn.list[index]<-sapply(cdn.list[index],rev)
	fram<-as.integer(sapply(cdn.list,function(x)x[1]))
	gene<-subseq(gene,fram+1)
	
	gene
}

ref<-readDNAStringSet(f.ref)
gff<-read.table(f.gff,sep="\t",stringsAsFactors=F,quote="")

gene<-gff2fa(gff,ref)
writeXStringSet(gene,f.out.gene)

prot<-translate(gene,if.fuzzy.codon="X")
i<-grep("\\*$",prot)
if(length(i)>0)prot[i]<-subseq(prot[i],1,width(prot[i])-1)
writeXStringSet(prot,f.out.prot)
