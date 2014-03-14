package de.helmholtz_muenchen.ibis.ngs.filtersummary.output;

import de.helmholtz_muenchen.ibis.ngs.filtersummary.genes.GeneList;
import de.helmholtz_muenchen.ibis.utils.lofs.Writer_Output;
import de.helmholtz_muenchen.ibis.ngs.filtersummary.vcf.LoFWorkflowResult;

public class SummaryFile {
	
	public static void summaryFile (LoFWorkflowResult r, String file){
		
		Writer_Output ws = new Writer_Output(file);
		
		//file counts
		String filecount= r.pindeldel.countsummary()+"\n"+
				r.pindelins.countsummary()+"\n"+
				r.gatksnps.countsummary()+"\n"+
				r.gatkindels.countsummary();
		
		System.out.println(filecount);
		ws.writeFileln(filecount);
		
		// all counts
		String all = getALLStatistics(r);
		System.out.println(all);
		ws.writeFileln(all);
		
		// lof counts
		String lofs = getLOFStatistics(r);
		System.out.println(lofs);
		ws.writeFileln(lofs);
		
		// snp-indel counts
		String snpindel = getIndelsnpCounts(r);
		System.out.println(snpindel);
		ws.writeFileln(snpindel);
		
		// gene counts
		System.out.println(getGeneCounts(r.lofgenes));
		ws.writeFileln(getGeneCounts(r.lofgenes));
		
		//filter counts
		String filter=getFilterStats(r);
		System.out.println(filter);
		ws.writeFileln(filter);
		
		//filtered genes
		System.out.println(getGeneCounts(r.lofgenes_filt));
		ws.writeFileln(getGeneCounts(r.lofgenes_filt));
		
		ws.closew();		
		
	}
	
	public static String getALLStatistics (LoFWorkflowResult r){
		String res ="All variants\n";
		res+="All variants\t"+r.varsall+"\n";
		res+="Duplicates\t"+r.duplall+"\n";
		res+="GATK\t"+r.gatkvar+"\n";
		res+="Pindel\t"+r.pindelvar+"\n";
		res+="DBSNP entries\t"+r.dbentries+"\n";
		res+="Homozygous\t"+r.homoall+"\n";
		return res;
	}
	
	public static String getLOFStatistics(LoFWorkflowResult r){
		String res = "LOF counts:\n";
		res+="Lines\t"+r.vars+"\n";
		res+="LoFs\t"+r.lofs+"\n";
		res+="Duplicates\t"+r.countdupl+"\n";
		res+="GATK\t"+r.gatklofs+"\n";
		res+="Pindel\t"+r.pindellofs+"\n";
		res+="DBSNP entries\t"+r.dblofentries+"\n";
		res+="Homozygous\t"+r.homo+"\n";
		res+="InsertionFS\t"+r.inFS+"\n";
		res+="DeletionFS\t"+r.delFS+"\n";
		res+="RemovedStop\t"+r.rstp+"\n";
		res+="PrematureStop\t"+r.pmstp+"\n";
		res+="SpliceOverlap\t"+r.so+"\n";
		res+="MultiExon\t"+r.multEx+"\n";
		return res;
	}
	
	public static String getIndelsnpCounts(LoFWorkflowResult r){
		String res="Snps\t"+r.snps+"\n";
		res+="No LoF Snps\t"+r.snpnolof+"\n";
		res+="Premature Stop Snp\t"+r.pmstp+"\n";
		res+="Removed Stop Snp\t"+r.rstp+"\n";
		res+="Splice Overlap Snp\t"+r.sosnp+"\n";
		res+="\n";
		res+="Insertions\t"+r.insertions+"\n";
		res+="No LoF Insertion\t"+r.insnolof+"\n";
		res+="Frame Shift Insertion\t"+r.inFS+"\n";
		res+="Splice Overlap Insertion\t"+r.soins+"\n";
		res+="\n";
		res+="Deletions\t"+r.deletions+"\n";
		res+="No LoF Deletion\t"+r.delnolof+"\n";
		res+="Frame Shift Deletion\t"+r.delFS+"\n";
		res+="Splice Overlap Deletion\t"+r.sodel+"\n";
		return res;
	}
	
	public static String getFilterStats(LoFWorkflowResult r){
		String res="Variants after filtering\n";
		res+=r.cf+"\n";
		res+="lof variants\t"+r.fvars+"\n";
		res+="gatk+pindel lofs\t"+r.fcountdupl+"\n";
		res+="gatk lofs\t"+r.fgatklofs+"\n";
		res+="pindel lofs\t"+r.fpindellofs+"\n";
		res+="homozygous lofs\t"+r.fhomo+"\n";
		res+="dbsnp lofs\t"+r.fdblofentries+"\n";
		res+="\n";
		res+="all lofs\t"+r.flofs+"\n";
		res+="insertionFS\t"+r.finFS+"\n";
		res+="deletionFS\t"+r.fdelFS+"\n";
		res+="removedStop\t"+r.frstp+"\n";
		res+="prematureStop\t"+r.frstp+"\n";
		res+="SpliceOverlap\t"+r.fso+"\n";
		res+="MultiExon\t"+r.fmultEx+"\n";
		return res;
	}
	
	public static String getGeneCounts(GeneList l){
		String res="LoF Genes\t"+l.genes+"\n";
		res+="GATK LoF Genes\t"+l.gatk+"\n";
		res+="Pindel LoF Genes\t"+l.pindel+"\n";
		res+="Full LoF Genes\t"+l.full+"\n";
		res+="Genes with OMIM entries\t"+l.omim+"\n";
		return res;
	}

}
