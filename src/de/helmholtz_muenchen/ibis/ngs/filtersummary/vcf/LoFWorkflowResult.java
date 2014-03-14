package de.helmholtz_muenchen.ibis.ngs.filtersummary.vcf;

import de.helmholtz_muenchen.ibis.ngs.filtersummary.filter.CovFilter;
import de.helmholtz_muenchen.ibis.ngs.filtersummary.genes.GeneInfo;
import de.helmholtz_muenchen.ibis.ngs.filtersummary.genes.GeneList;
import de.helmholtz_muenchen.ibis.ngs.filtersummary.interactions.ReactomePPI;

import java.util.Iterator;
import java.util.LinkedList;

import de.helmholtz_muenchen.ibis.ngs.filtersummary.output.FilesGenes;
import de.helmholtz_muenchen.ibis.ngs.filtersummary.output.PlotFilesVar;
import de.helmholtz_muenchen.ibis.ngs.filtersummary.output.SummaryFile;
import de.helmholtz_muenchen.ibis.ngs.filtersummary.output.VariantFile;
import de.helmholtz_muenchen.ibis.ngs.filtersummary.pathways.WikiPathway;

public class LoFWorkflowResult {
	
	//4 vcf files, output of LoF workflow
	public VCFFilePindelDel pindeldel;
	public VCFFilePindelIns pindelins;
	public VCFFileGatkSnp gatksnps;
	public VCFFileGatkIndel gatkindels;
	
	private String outbase;
	
	public LinkedList<VCFLine> varlist;
	
	public LinkedList<VCFLine> loflist;
	
	public LinkedList<VCFLine> loflist_filt;
	
	public CovFilter cf;
	
	public int varsall;		//vcf line all variants
	public int duplall;		// duplcates all variants
	public int homoall;		// homozygous all variants
	public int pindelvar;		//variants found by pindel
	public int gatkvar;		//varaints found by gatk
	public int dbentries;
	
	public int vars;		//vcf lines lof
	public int countdupl;	// varaints found by gatk and pindel -> count lines
	public int gatklofs;
	public int pindellofs;
	public int homo;		// homozygous variants -> count lines
	public int dblofentries;
	
	public int lofs;		// vat lofs annotations, more than one annotation per line possible
	public int inFS=0;
	public int delFS=0;
	public int multEx=0;
	public int so=0;
	public int rstp =0;
	public int pmstp =0;

	public int insertions=0;	// analyse proportion of insertions/deletions/snps that are some lofs
	public int deletions=0;
	public int snps=0;
	public int soins=0;
	public int sodel=0;
	public int sosnp=0;
	public int snpnolof=0;
	public int delnolof=0;
	public int insnolof=0;
	
	// counts for filtered lofs
	public int fvars;
	public int fcountdupl;
	public int fgatklofs;
	public int fpindellofs;
	public int fhomo;
	public int fdblofentries;
	
	public int flofs;
	public int finFS;
	public int fdelFS;
	public int fmultEx;
	public int fso;
	public int frstp;
	public int fpmstp;
	
	
	public GeneList lofgenes;
	public GeneList lofgenes_filt;
	
	public GeneInfo geneinfo;
	public WikiPathway wikipathway;
	public ReactomePPI reactomeppi;
	
	public String ppifile;
	
	
	public LoFWorkflowResult (String pindeldel, String pindelins, String gatksnps, String gatkindels, String outbase, String geneinfofile, String pathwayfile, String ppifile){
		this.pindeldel=new VCFFilePindelDel(pindeldel);
		this.pindelins=new VCFFilePindelIns(pindelins);
		this.gatksnps=new VCFFileGatkSnp(gatksnps);
		this.gatkindels=new VCFFileGatkIndel(gatkindels);
		
		this.outbase=outbase;
		
		loflist = new LinkedList<VCFLine>();
		varlist = new LinkedList<VCFLine>();
		
		if(geneinfofile!=null){
			geneinfo= new GeneInfo(geneinfofile);
		}
		else{
			geneinfo=null;
		}
		
		if(pathwayfile!=null){
			wikipathway= new WikiPathway(pathwayfile);
		}
		else{
			wikipathway=null;
		}
		
		if(ppifile!=null){
			reactomeppi=new ReactomePPI(ppifile);
		}
		else{
			reactomeppi=null;
		}
		
		cf = new CovFilter();
	}
	
	// general count statistics
	
	public void load() throws Exception{
		pindeldel.load();
		pindelins.load();
		gatksnps.load();
		gatkindels.load();
		
		pindelvar=pindeldel.getVarCount()+pindelins.getVarCount();
		gatkvar=gatksnps.getVarCount()+gatkindels.getVarCount();
		dbentries=gatksnps.getdbsnpcount()+gatkindels.getdbsnpcount();
	}
	
	// merge vcffiles 
	
	public void mergeVars (){
		
		gatkindels.VarIterator();
		gatksnps.VarIterator();
		pindeldel.VarIterator();
		pindelins.VarIterator();
		
		VCFLine indel = gatkindels.nextVariant();
		VCFLine snp = gatksnps.nextVariant();
		VCFLine deletion = pindeldel.nextVariant();
		VCFLine insertion = pindelins.nextVariant();
		
		while (!(indel==null && snp==null && deletion==null && insertion==null) ){
			
			if(indel!=null){
				if(indel.equals(snp)){
					duplall++;
					indel = gatkindels.nextVariant();
					System.out.println("gatk indel==gatk snp");
					continue;
				}
				
				if(indel.equals(deletion)){
					duplall++;
					deletion = pindeldel.nextVariant();
					indel.setTool(VCFFile.GATK+VCFFile.PINDEL);
					continue;
				}
				
				if(indel.equals(insertion)){
					duplall++;
					insertion = pindelins.nextVariant();
					indel.setTool(VCFFile.GATK+VCFFile.PINDEL);
					continue;
				}
			}
			
			if(snp != null){
				if(snp.equals(deletion)){
					duplall++;
					deletion = pindeldel.nextVariant();
					snp.setTool(VCFFile.GATK+VCFFile.PINDEL);
					System.out.println("snp==deletion");
					continue;
				}
				
				if(snp.equals(insertion)){
					duplall++;
					insertion = pindelins.nextVariant();
					snp.setTool(VCFFile.GATK+VCFFile.PINDEL);
					System.out.println("snp==insertion");
					continue;
				}
			}
			
			if(deletion != null){
				if(deletion.equals(insertion)){
					duplall++;
					deletion = pindeldel.nextVariant();
					System.out.println("insertion==deletion");
					continue;
				}
			}
			
			
			VCFLine min = VCFFile.minimum (indel, snp, deletion, insertion);
			varlist.add(min);
			//System.out.println(min.toString());
			
			varsall++;
			
			if(min.getSample().startsWith("1/1") || min.getSample().startsWith("1|1")){
				homoall++;
			}
			
			if(min==indel){
				indel=gatkindels.nextVariant();
			}
			
			if(min==snp){
				snp = gatksnps.nextVariant();
			}
			
			if(min==deletion){
				deletion = pindeldel.nextVariant();
			}
			
			if(min==insertion){
				insertion = pindelins.nextVariant();
			}
			
		}
		
	}
	
	public void mergeLoFs (){
		
		gatkindels.LoFIterator();
		gatksnps.LoFIterator();
		pindeldel.LoFIterator();
		pindelins.LoFIterator();
		
		VCFLine indel = gatkindels.nextLoFVariant();
		VCFLine snp = gatksnps.nextLoFVariant();
		VCFLine deletion = pindeldel.nextLoFVariant();
		VCFLine insertion = pindelins.nextLoFVariant();
		
		while (!(indel==null && snp==null && deletion==null && insertion==null) ){
			
			if(indel!=null){
				if(indel.equals(snp)){
					countdupl++;
					indel = gatkindels.nextLoFVariant();
					System.out.println("gatk indel==gatk snp");
					continue;
				}
				
				if(indel.equals(deletion)){
					countdupl++;
					deletion = pindeldel.nextLoFVariant();
					indel.setTool(VCFFile.GATK+VCFFile.PINDEL);
					continue;
				}
				
				if(indel.equals(insertion)){
					countdupl++;
					insertion = pindelins.nextLoFVariant();
					indel.setTool(VCFFile.GATK+VCFFile.PINDEL);
					continue;
				}
			}
			
			if(snp != null){
				if(snp.equals(deletion)){
					countdupl++;
					deletion = pindeldel.nextLoFVariant();
					snp.setTool(VCFFile.GATK+VCFFile.PINDEL);
					System.out.println("snp==deletion");
					continue;
				}
				
				if(snp.equals(insertion)){
					countdupl++;
					insertion = pindelins.nextLoFVariant();
					snp.setTool(VCFFile.GATK+VCFFile.PINDEL);
					System.out.println("snp==insertion");
					continue;
				}
			}
			
			if(deletion != null){
				if(deletion.equals(insertion)){
					countdupl++;
					deletion = pindeldel.nextLoFVariant();
					System.out.println("insertion==deletion");
					continue;
				}
			}
			
			
			VCFLine min = VCFFile.minimum (indel, snp, deletion, insertion);
			loflist.add(min);
			//System.out.println(min.toString());
			
			vars++;
			
			if(min.getTool()==VCFFile.GATK){
				gatklofs++;
			}
			else if (min.getTool()==VCFFile.PINDEL){
				pindellofs++;
			}
			else{
				gatklofs++;
				pindellofs++;
			}
			
			if(!min.getID().equals(".")){
				dblofentries++;
			}
			
			if(min.getSample().startsWith("1/1") || min.getSample().startsWith("1|1")){
				homo++;
			}
			if(min.getspliceOverlap()){
				so++;
				lofs++;
			}
			if(min.getInsertionFS()){
				inFS++;
				lofs++;
			}
			if(min.getDeletionFS()){
				delFS++;
				lofs++;
			}
			if(min.getPrematureStop()){
				pmstp++;
				lofs++;
			}
			if(min.getRemovedStop()){
				rstp++;
				lofs++;
			}
			if(min.getMultiExon()){
				multEx++;
				lofs++;
			}

			
			
			if(min==indel){
				indel=gatkindels.nextLoFVariant();
			}
			
			if(min==snp){
				snp = gatksnps.nextLoFVariant();
			}
			
			if(min==deletion){
				deletion = pindeldel.nextLoFVariant();
			}
			
			if(min==insertion){
				insertion = pindelins.nextLoFVariant();
			}
			
		}
		
	}
	
	public void countindelsnp(){
		
		Iterator<VCFLine> it = varlist.iterator();
		VCFLine l;
		
		while(it.hasNext()){
			l=it.next();
			
			if(l.getSnpindel()==VCFLine.SNP){
				snps++;
				if(l.getInfo().matches(".*spliceOverlap.*")){
					sosnp++;
				}
				else if(!l.getInfo().matches(".*prematureStop.*") && !l.getInfo().matches(".*removedStop.*")){
					snpnolof++;
				}
			}
			else if(l.getSnpindel()==VCFLine.INS){
				insertions++;
				if(l.getInfo().matches(".*spliceOverlap.*")){
					soins++;
				}
				else if(!l.getInfo().matches(".*insertionFS.*")){
					insnolof++;
				}
			}
			else if(l.getSnpindel()==VCFLine.DEL){
				deletions++;
				if(l.getInfo().matches(".*spliceOverlap.*")){
					sodel++;
				}
				else if(!l.getInfo().matches(".*deletionFS.*")){
					delnolof++;
				}
			}
			else if(l.getSnpindel()==VCFLine.INDEL){
				insertions++;
				deletions++;
				if(l.getInfo().matches(".*spliceOverlap.*")){
					soins++;
					sodel++;
				}
				else if(!l.getInfo().matches(".*deletionFS.*")){
					delnolof++;
				}
				else if(!l.getInfo().matches(".*insertionFS.*")){
					insnolof++;
				}
			}
		}
	}
	
	
	public void loffiltering(boolean summary){
		loflist_filt = cf.filter(this.loflist);
		
		if(summary){
			fvars=loflist_filt.size();
			
			Iterator<VCFLine> i = loflist_filt.iterator();
			while(i.hasNext()){
				VCFLine v = i.next();
				
				//get tool counts
				if(v.getTool()==VCFFile.PINDEL+VCFFile.GATK){
					fcountdupl++;
					fgatklofs++;
					fpindellofs++;
				}
				else if(v.getTool()==VCFFile.PINDEL){
					fpindellofs++;
				}
				else if(v.getTool()==VCFFile.GATK){
					fgatklofs++;
				}
				
				//count homozygous variants
				if(v.getGenotype().equals("1/1") ||v.getGenotype().equals("1|1")){
					fhomo++;
				}
				
				//count dbsnp entries
				if(!v.getID().equals(".")){
					fdblofentries++;
				}
				
				//count lof types
				if(v.getDeletionFS()){
					flofs++;
					fdelFS++;
				}
				if(v.getInsertionFS()){
					flofs++;
					finFS++;
				}
				if(v.getMultiExon()){
					flofs++;
					fmultEx++;
				}
				if(v.getPrematureStop()){
					flofs++;
					fpmstp++;
				}
				if(v.getRemovedStop()){
					flofs++;
					frstp++;
				}
				if(v.getspliceOverlap()){
					flofs++;
					fso++;
				}
			}		
		}
	}
	
	public void plots(){
		
		//files for plots
		//statistics
		PlotFilesVar.writeCountTable(this, outbase+"counttable.txt");		// bar plot counting variants
		PlotFilesVar.writeTypeTable(this, outbase+"typetable.txt");			// pie plot types of lof
		PlotFilesVar.writeSnpIndelTable(this, outbase+"snpindelcount.txt");	// stacked bar plot proportion of lof insertions/deletions
		
		// dbsnp vs no dbsnp
		PlotFilesVar.writeQual(this, outbase+"nodbentry.txt", outbase+"dbentry.txt");		// ecdf for gatk scores and dbsnp entries
		PlotFilesVar.writeCovgatk(this, outbase+"dbcov.txt", outbase+"nodbcov.txt");		// ecdf for coverage and dbsnp entries
		
		// gatk vs pindel
		PlotFilesVar.writeLengthDel(this, outbase+"delgatk.txt", outbase+"delpindel.txt");	// length of deletions
		PlotFilesVar.writeLengthIns(this, outbase+"insgatk.txt", outbase+"inspindel.txt");	// length of insertions
		PlotFilesVar.writeQualgatkpindel(this, outbase+"indelgatkqual.txt", outbase+"indelgatkpindelqual.txt");	// quality scores of gatk indels
		PlotFilesVar.writeIndelCounts(this, outbase+"indeltooltable.txt");	// indel counts for pindel and gatk
		
		//genes
		FilesGenes.writeGeneList(lofgenes, outbase+"genes.txt");
		FilesGenes.writeLofsperGene(lofgenes, outbase+"genelistcount.txt");	// lofs per gene
		FilesGenes.writeGeneCountTable(lofgenes, outbase+"genecount.txt");	// gene counts that are also listed in summary file
		
		//interactions of lof genes compared to all genes in reactome data set
		FilesGenes.writeGeneDegree(lofgenes, reactomeppi, outbase+"ppideg.txt", outbase+"allppideg.txt", 0);	// all interactions
		FilesGenes.writeGeneDegree(lofgenes, reactomeppi, outbase+"cppideg.txt", outbase+"allcppideg.txt", 1);	// only direct complexes
		FilesGenes.writeGeneDegree(lofgenes, reactomeppi, outbase+"rppideg.txt", outbase+"allrppideg.txt", 2);	// only direct reaction
		
		//pathways from wikipathways in which lof genes are involved
		FilesGenes.writePathways(lofgenes, wikipathway, outbase+"pathways.txt");
		
	}
	
	public void runAnalysis(boolean summary, int tables, boolean plots) throws Exception{
		
		//variant analysis
		this.load();
		this.mergeVars();
		this.mergeLoFs();
		
		// gene analysis
		lofgenes=new GeneList(loflist);
		lofgenes.countGenes(geneinfo);	
		
		// filter
		loffiltering(summary);
		lofgenes_filt=new GeneList(loflist_filt);
		lofgenes_filt.countGenes(geneinfo);

		
		if(summary){
			this.countindelsnp();
			SummaryFile.summaryFile(this,outbase+"summary.txt");
		}
		
		String vcfs = gatksnps.getPath()+"\t"+gatkindels.getPath()+"\t"+pindeldel.getPath()+"\t"+pindelins.getPath();
		if(tables!=0){
			VariantFile.createVarFile(vcfs, varlist, geneinfo, reactomeppi, wikipathway, outbase+"all.csv", false, tables);
			VariantFile.createVarFile(vcfs, loflist, geneinfo, reactomeppi, wikipathway, outbase+"lofs.csv", true, tables);
			VariantFile.createVarFile(vcfs, loflist_filt, geneinfo, reactomeppi, wikipathway, outbase+"filteredlofs.csv", true, tables);
		}
		
		if(plots){
			this.plots();
		}
	}

}
