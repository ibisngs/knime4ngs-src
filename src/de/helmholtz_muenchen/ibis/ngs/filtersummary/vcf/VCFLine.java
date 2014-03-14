package de.helmholtz_muenchen.ibis.ngs.filtersummary.vcf;

import de.helmholtz_muenchen.ibis.ngs.filtersummary.genes.Gene;

import java.util.Iterator;
import java.util.LinkedList;

public class VCFLine {
	
	public static final int SNP=1;
	public static final int INS=2;
	public static final int DEL=3;
	public static final int INDEL=4;
	
	private String chrom; 		//0
	private int pos;			//1
	private String id;			//2
	private String ref;			//3
	private String alt;			//4
	private String qual;		//5
	//String filter;
	private String info;		//7
	//String format;
	private String sample;		//9
	
	// retrieve from sample field
	private String genotype;
	private int coverage;
	private int varreads;
	
	private int tool;			// gatk, pindel or gatk+pindel
	private int snpindel;		// snp, insertion, deletion, indel or nothing
	
	// lof types
	private boolean spliceOverlap;
	private boolean insertionFS;
	private boolean deletionFS;
	private boolean prematureStop;
	private boolean removedStop;
	private boolean multiExon;
	// no lof types
	private boolean syn;
	private boolean nonsyn;
	private boolean insertionNFS;
	private boolean deletionNFS;
	private boolean startOverlap;
	private boolean endOverlap;
	
	private LinkedList<Gene> affgenes;
	
	
	public VCFLine (String line, int tool, int snpindel) throws Exception{
		
		
		String [] split= line.split("\t");
		
		if(split.length!=10){
			throw new Exception("Wrong fromat!"+line);
		}
		
		chrom=split[0];
		pos=Integer.parseInt(split[1]);
		id=split[2];
		ref=split[3];
		alt=split[4];
		qual=split[5];
		info=split[7];
		sample=split[9];
		
		this.tool=tool;
		this.snpindel=snpindel;
		

		String [] sa = sample.split(":");
		genotype=sa[0];
		String [] s = sa[1].split(",");
		varreads=Integer.parseInt(s[1]);
		int refreads=Integer.parseInt(s[0]);
		coverage=varreads+refreads;
		
		this.setMutationType();
		
		
		// create list with genes
		affgenes= new LinkedList<Gene>();
		
		this.findGenes();
		
	}
	
	private void setMutationType(){
		
		if(info.matches(".*spliceOverlap.*")){
			spliceOverlap=true;
		}
		
		if(info.matches(".*insertionFS.*")){
			insertionFS=true;
		}
		
		if(info.matches(".*deletionFS.*")){
			deletionFS=true;
		}
		
		if(info.matches(".*prematureStop.*")){
			prematureStop=true;
		}
		
		if(info.matches(".*removedStop.*")){
			removedStop=true;
		}
		
		if(info.matches(".*multiExon.*")){
			multiExon=true;
		}

		if(info.matches(".*synonymous.*")){
			if(info.matches(".*nonsynonymous.*")){
				nonsyn=true;
			}
			else{
				syn=true;
			}
		}
		
		if(info.matches(".*insertionNFS.*")){
			insertionNFS=true;
		}
		
		if(info.matches(".*deletionNFS.*")){
			deletionNFS=true;
		}
		
		if(info.matches(".*startOverlap.*")){
			startOverlap=true;
		}
		if(info.matches(".*endOverlap.*")){
			endOverlap=true;
		}
	}
	
	public void findGenes() throws Exception{
		
		String splitinfo [] = info.split("VA=");
		
		if(splitinfo.length!=2){
			throw new Exception("Missing VAT annotation: "+info);
		}
		
		String splitvat [] = splitinfo[1].split(",");
		for(int i=0; i<splitvat.length;i++){
			Gene g = new Gene (splitvat[i]);
			affgenes.add(g);
		}
		
	}
	
	public String toString(){
		
		String caller="";
		
		if(tool==VCFFile.PINDEL){
			caller="pindel";
		}
		else if (tool ==VCFFile.GATK){
			caller="gatk";
		}
		else{
			caller="gatk+pindel";
		}
		
		String vartype="unknown type";
		if(snpindel==VCFLine.DEL){
			vartype="deletion";
		}
		else if(snpindel==VCFLine.INS){
			vartype="insertion";
		}
		else if(snpindel==VCFLine.INDEL){
			vartype="indel";
		}
		else if(snpindel==VCFLine.SNP){
			vartype="snp";
		}
		
		return caller+"\t"+vartype+"\t"+chrom+"\t"+pos+"\t"+id+"\t"+ref+"\t"+alt+"\t"+qual+"\t"+info+"\t"+sample;
		
	}

	
	public void setTool(int toolnum){
		
		if(toolnum==VCFFile.GATK+VCFFile.PINDEL){
			this.tool=VCFFile.GATK+VCFFile.PINDEL;
		}
		
	}
	
	public void setSnpindel(int snpindel){
		if(snpindel==SNP || snpindel==DEL || snpindel==INS || snpindel==INDEL){
			this.snpindel=snpindel;
		}
	}
	

	
	public static VCFLine min (VCFLine a, VCFLine b){
		
		VCFLine min =null;
		
		if(a == null && b==null){
			return min;
		}
		
		if(a == null){
			return b;
		}
		
		if(b==null){
			return a;
		}
		
		int comp = new Chromosome(a.chrom).compareTo(new Chromosome(b.chrom));
		if(comp<0){
			min =a;
		}
		if(comp > 0){
			min = b;
		}
		if(comp == 0){
			
			if(a.pos< b.pos){
				min =a;
			}
			if(a.pos > b.pos){
				min =b;
			}
			if(a.pos == b.pos){
				min =a;
			}
		}
		
		return min;
		
	}

	public boolean equals(VCFLine x){
		
		if(x==null){
			return false;
		}
		
		int comp = new Chromosome(this.chrom).compareTo(new Chromosome(x.chrom));
		
		if(comp ==0 && this.pos==x.pos && this.ref.equals(x.getRef()) && this.alt.equals(x.alt)){
			return true;
		}
		
		return false;
		
	}
	
	
	public String getChrom(){
		return chrom;
	}
	
	public int getPosition(){
		return pos;
	}
	
	public int getLength(){
		if(this.snpindel==SNP){
			return 1;
		}
		return Math.abs(ref.length()-alt.length());
	}
	
	public String getSample(){
		return sample;
	}
	
	public String getInfo(){
		return info;
	}
	
	public String getRef(){
		return ref;
	}
	
	public String getAlt(){
		return alt;
	}
	
	public int getTool(){
		return tool;
	}
	
	public String getToolName(){
		String tool="";
		if(this.tool==VCFFile.GATK){
			tool="GATK";
		}
		else if(this.tool==VCFFile.PINDEL){
			tool="Pindel";
		}
		else if(this.tool==VCFFile.GATK+VCFFile.PINDEL){
			tool="GATK+Pindel";
		}
		return tool;
	}
	
	public String getID(){
		return id;
	}
	
	public String getQual(){
		return qual;
	}
	
	public int getSnpindel(){
		return snpindel;
	}
	
	public String getVariantType(){
		String type="";
		if(this.snpindel==SNP){
			type="SNP";
		}
		else if(this.snpindel==DEL){
			type="Deletion";
		}
		else if(this.snpindel==INS){
			type="Insertion";
		}
		else if(this.snpindel==INDEL){
			type="Indel";
		}
		return type;
	}
	
	public int getCov(){
		return coverage;
	}
	
	public int getSuppReads(){
		return varreads;
	}
	
	public String getGenotype(){
		return genotype;
	}
	
	public LinkedList<Gene> getGenes(){
		return this.affgenes;
	}
	
	public String getGenes (boolean lof, int info){
		String geneids="";
		Iterator<Gene> i = affgenes.iterator();
		while(i.hasNext()){
			Gene g = i.next();
			if(lof){
				if(!g.getType().equals("prematureStop") && !g.getType().equals("spliceOverlap") && !g.getType().equals("removedStop") && 
					!g.getType().equals("deletionFS") && !g.getType().equals("insertionFS") && !g.getType().equals("multiExonHit")){
					continue;
				}	
			}
			if(!geneids.equals("")){
				geneids+=",";
			}
			if(info==0){
				geneids+=g.getIDnoExt();
			}
			else if(info==1){
				geneids+=g.getName();
			}
			else if(info==2){
				geneids+=g.getAffTranscripts();
			}
			else if(info==3){
				geneids+=g.getTranscripts();
			}
		}
		return geneids;
	}
	
	public boolean getspliceOverlap(){
		return this.spliceOverlap;
	}
	
	public boolean getInsertionFS(){
		return this.insertionFS;
	}
	
	public boolean getDeletionFS(){
		return this.deletionFS;
	}
	
	public boolean getPrematureStop(){
		return this.prematureStop;
	}
	
	public boolean getRemovedStop(){
		return this.removedStop;
	}
	
	public boolean getMultiExon(){
		return this.multiExon;
	}
	
	public boolean getSyn(){
		return this.syn;
	}
	
	public boolean getNonsyn(){
		return this.nonsyn;
	}
	
	public boolean getInsertionNFS(){
		return this.insertionNFS;
	}
	
	public boolean getDeletionNFS(){
		return this.deletionNFS;
	}
	
	public boolean getStartoverlap(){
		return this.startOverlap;
	}
	
	public boolean getEndoverlap(){
		return this.endOverlap;
	}
	
	public String getLofs( boolean lof){
		String loftype="";
		if(this.insertionFS){
			if(!loftype.equals("")){
				loftype+=",";
			}
			loftype+="Insertion Frameshift";
		}
		if(this.deletionFS){
			if(!loftype.equals("")){
				loftype+=",";
			}
			loftype+="Deletion Frameshift";
		}
		if(this.removedStop){
			if(!loftype.equals("")){
				loftype+=",";
			}
			loftype+="Removed Stop";
		}
		if(this.prematureStop){
			if(!loftype.equals("")){
				loftype+=",";
			}
			loftype+="Premature Stop";
		}
		if(this.spliceOverlap){
			if(!loftype.equals("")){
				loftype+=",";
			}
			loftype+="Splice Overlap";
		}
		if(this.multiExon){
			if(!loftype.equals("")){
				loftype+=",";
			}
			loftype+="Multi-Exon Hit";
		}
		if(!lof && this.deletionNFS){
			if(!loftype.equals("")){
				loftype+=",";
			}
			loftype+="Deletion without frameshift";
		}
		if(!lof && this.insertionNFS){
			if(!loftype.equals("")){
				loftype+=",";
			}
			loftype+="Insertion without frameshift";
		}
		if(!lof && this.startOverlap){
			if(!loftype.equals("")){
				loftype+=",";
			}
			loftype+="Start Overlap";
		}
		if(!lof && this.endOverlap){
			if(!loftype.equals("")){
				loftype+=",";
			}
			loftype+="End Overlap";
		}
		if(!lof && this.syn){
			if(!loftype.equals("")){
				loftype+=",";
			}
			loftype+="Synonymous";
		}
		if(!lof && this.nonsyn){
			if(!loftype.equals("")){
				loftype+=",";
			}
			loftype+="Non-synonymous";
		}
		return loftype;
	}
}
