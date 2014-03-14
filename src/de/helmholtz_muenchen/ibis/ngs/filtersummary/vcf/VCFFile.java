package de.helmholtz_muenchen.ibis.ngs.filtersummary.vcf;

import java.util.Iterator;
import java.util.LinkedList;

import de.helmholtz_muenchen.ibis.utils.lofs.FileInputReader;

public class VCFFile {
	
	public static final int GATK = 1;
	public static final int PINDEL = 2;
	
	private String path;
	private FileInputReader fir;
	
	protected LinkedList<VCFLine> entries = new LinkedList<VCFLine>();
	protected LinkedList<VCFLine> LoFs = new LinkedList<VCFLine>();
	
	private Iterator<VCFLine> itentries = null;
	private Iterator<VCFLine> itlofs = null;
	
	private int tool;
	
	private int variants =0;
	private int lofs =0;
	
	//lofs
	private int spliceOverlap =0;
	private int insertionFS =0;
	private int deletionFS =0;
	private int prematureStop =0;
	private int removedStop =0;
	private int multiExon =0;
	
	private int syn=0;
	private int nonsyn=0;
	private int insertionNFS=0;
	private int deletionNFS=0;
	private int startOverlap=0;
	private int endOverlap=0;
	
	public VCFFile(String path, int tool){
		this.path = path;
		
		try{
		fir = new FileInputReader(path);
		}
		catch(Exception e){
			e.printStackTrace();
		}
		
		if(tool!=GATK && tool!=PINDEL){
			System.out.println("unkown tool");
		}
		if(tool==GATK){
			this.tool =GATK;
		}
		if(tool==PINDEL){
			this.tool = PINDEL;
		}
	}
	
	public int getType(VCFLine line){
		return 0;
	}
	
	public void load () throws Exception{
		String line;
		
		while((line=fir.read())!=null){

			
			if(!line.startsWith("#")){
				
				//remove calls from pindel with genotype 0-0
				String [] split = line.split("\t");
				if(!split[split.length-1].startsWith("0/0") && !split[split.length-1].startsWith("0|0")){
					

					variants++;
					
					VCFLine finalvar = new VCFLine(line, tool, 0);
					
					//decide if variant is snp, insertion or deletion
					int type=getType(finalvar);
					if(type==0){
						System.out.println("Cannot identify variant type");
						System.out.println(line);
					}					
					finalvar.setSnpindel(type);
					
					entries.add(finalvar);
										
					//some variants annotated with 2 or more lof types -> only once in LoFs list, but count as 2 or more lofs
					boolean lof=false;
					
					if(finalvar.getspliceOverlap()){
						spliceOverlap++;
						lofs++;
						lof=true;
					}
					
					if(finalvar.getInsertionFS()){
						lofs++;
						insertionFS++;
						lof=true;
					}
					
					if(finalvar.getDeletionFS()){
						lofs++;
						deletionFS++;
						lof=true;
					}
					
					if(finalvar.getPrematureStop()){
						lofs++;
						prematureStop++;
						lof=true;
					}
					
					if(finalvar.getRemovedStop()){
						lofs++;
						removedStop++;
						lof=true;
					}
					
					if(finalvar.getMultiExon()){
						lofs++;
						multiExon++;
						lof=true;
					}

					// add to lof list
					if(lof){
						LoFs.add(finalvar);
					}
					
					//count non-lof types
					
					if(finalvar.getSyn()){
						syn++;
					}
					
					if(finalvar.getNonsyn()){
						nonsyn++;
					}
					
					if(finalvar.getInsertionNFS()){
						insertionNFS++;
					}
					
					if(finalvar.getDeletionNFS()){
						deletionNFS++;
					}
					
					if(finalvar.getStartoverlap()){
						startOverlap++;
					}
					if(finalvar.getEndoverlap()){
						endOverlap++;
					}

				}
			}
		}
	}
	
	// get summary of vcf file data
	
	public int getVarCount(){
		return this.variants;
	}
	
	public String countsummary (){
		
		String res=path+"\n";
		res+="Variants\t"+variants+"\n";
		res+="LoFs\t"+lofs+"\n";
		res+="SpliceOverlap\t"+spliceOverlap+"\n";
		res+="InsertionFS\t"+ insertionFS+"\n";
		res+="DeletionFS\t"+deletionFS+"\n";
		res+="PrematureStop\t"+prematureStop+"\n";
		res+="RemovedStop\t"+removedStop+"\n";
		res+="MultiExon\t"+multiExon+"\n";
		res+="Synonymous\t"+syn+"\n";
		res+="Nonsynonymous\t"+nonsyn+"\n";
		res+="InsertionNFS\t"+insertionNFS+"\n";
		res+="DeletionNFS\t"+deletionNFS+"\n";
		res+="StartOverlap\t"+startOverlap+"\n";
		res+="EndOverlap\t"+endOverlap+"\n";
		
		return res;
	}
	
	// methods for merging vcf file data
	
	public VCFLine nextVariant(){
		
		if(itentries.hasNext()){
			return itentries.next();
		}
		else{
			return null;
		}
	}
	
	public void VarIterator(){
		itentries =entries.iterator();
	}
	
	public VCFLine nextLoFVariant(){

		if(itlofs.hasNext()){
			return itlofs.next();
		}
		else{
			return null;
		}
	}
	
	public void LoFIterator(){
		itlofs = LoFs.iterator();
	}
	
	public static VCFLine minimum (VCFLine a, VCFLine b, VCFLine c, VCFLine d){

		return VCFLine.min(VCFLine.min(a, b), VCFLine.min(c, d));
		
	}	
	
	public String getPath(){
		return this.path;
	}

}
