package de.helmholtz_muenchen.ibis.ngs.filtersummary.output;

import java.util.Iterator;

import de.helmholtz_muenchen.ibis.utils.lofs.Writer_Output;
import de.helmholtz_muenchen.ibis.ngs.filtersummary.vcf.LoFWorkflowResult;
import de.helmholtz_muenchen.ibis.ngs.filtersummary.vcf.VCFFile;
import de.helmholtz_muenchen.ibis.ngs.filtersummary.vcf.VCFLine;

public class PlotFilesVar {
	
	//bar chart variant counts
		public static void writeCountTable(LoFWorkflowResult r, String path){
			
			Writer_Output wo = new Writer_Output(path);
			
			wo.writeFileln("Variants\tGATK\tPindel\tDBSNPentries\tHomozygous"+
					"\tLoFs\tGATKLoFs\tPindelLoFs\tDBSNPentriesLoFs\tHomozygousLoFs");
			
			wo.writeFileln(r.varsall+"\t"+r.gatkvar+"\t"+r.pindelvar+"\t"+r.dbentries+"\t"+r.homoall+
					"\t"+r.vars+"\t"+r.gatklofs+"\t"+r.pindellofs+"\t"+r.dblofentries+"\t"+r.homo);
			
			wo.closew();
			
		}
		
	//pie chart lof types
		public static void writeTypeTable(LoFWorkflowResult r,String path){
			
			Writer_Output wo = new Writer_Output(path);
			
			wo.writeFileln("InsertionFS\tDeletionFS\tRemovedStop\tPrematureStop\tSpliceOverlap\tMultiExon");
			wo.writeFileln(r.inFS+"\t"+r.delFS+"\t"+r.rstp+"\t"+r.pmstp+"\t"+r.so+"\t"+r.multEx);
			
			wo.closew();
			
		}
		
	//ecdf of qualscore dbsnp vs. no dbsnp
		public static void writeQual(LoFWorkflowResult r, String pathnoentry, String pathentry){
			
			Writer_Output w1 = new Writer_Output(pathnoentry);
			Writer_Output w2 = new Writer_Output(pathentry);
			
			Iterator<VCFLine> it = r.varlist.iterator();
			
			VCFLine l;
			while(it.hasNext()){
				
				l=it.next();
				
				if(l.getTool()==VCFFile.GATK){
					if(l.getID().equals(".")){
						w1.writeFileln(l.getQual());
					}
					else{
						w2.writeFileln(l.getQual());
					}
				}
			}
			
			w1.closew();
			w2.closew();
			
		}
		
	// histogram of deletion length pindel vs. gatk
		public static void writeLengthDel(LoFWorkflowResult r, String pathgatk, String pathpindel){
			
			Writer_Output w1 = new Writer_Output(pathgatk);
			Writer_Output w2 = new Writer_Output(pathpindel);
			
			Iterator<VCFLine> it = r.varlist.iterator();
			
			VCFLine l;
			
			while(it.hasNext()){
				
				l=it.next();
				
				int reflen=l.getRef().length();
				int altlen=l.getAlt().length();
				
				int length=reflen-altlen;
				
				//select deletions
				if(!(reflen==1 && altlen==1) && length>0){
					if(l.getTool()==VCFFile.GATK){
						w1.writeFileln(Integer.toString(length));
					}
					else if(l.getTool()==VCFFile.PINDEL){
						w2.writeFileln(Integer.toString(length));
					}
					else{
						w1.writeFileln(Integer.toString(length));
						w2.writeFileln(Integer.toString(length));
					}
				}

			}
			
			w1.closew();
			w2.closew();
			
		}
		
	// histogram of insertion length pindel vs. gatk
		public static void writeLengthIns(LoFWorkflowResult r, String pathgatk, String pathpindel){
			
			Writer_Output w1 = new Writer_Output(pathgatk);
			Writer_Output w2 = new Writer_Output(pathpindel);
			
			Iterator<VCFLine> it = r.varlist.iterator();
			
			VCFLine l;
			
			while(it.hasNext()){
				
				l=it.next();
				
				int reflen=l.getRef().length();
				int altlen=l.getAlt().length();
				
				int length=altlen-reflen;
				
				//select insertions
				if(!(reflen==1 && altlen==1) && length>0){
					if(l.getTool()==VCFFile.GATK){
						w1.writeFileln(Integer.toString(length));
					}
					else if(l.getTool()==VCFFile.PINDEL){
						w2.writeFileln(Integer.toString(length));
					}
					else{
						w1.writeFileln(Integer.toString(length));
						w2.writeFileln(Integer.toString(length));
					}
				}

			}
			
			w1.closew();
			w2.closew();
		}
		
	// bar plot lof proportion of insertions/ deletions
		public static void writeSnpIndelTable(LoFWorkflowResult r, String file){
			Writer_Output w = new Writer_Output(file);
			w.writeFileln("#NolofSnp\tprestopSnp\tremstopSnp\tsoSnp\tNolofIns\tFsIns\tsoIns\tNolosDel\tFSDel\tsoDel");
			w.writeFileln(Integer.toString(r.snpnolof)+"\t"+Integer.toString(r.pmstp)+"\t"+Integer.toString(r.rstp)+"\t"+Integer.toString(r.sosnp)
					+"\t"+Integer.toString(r.insnolof)+"\t"+Integer.toString(r.inFS)+"\t"+Integer.toString(r.soins)
					+"\t"+Integer.toString(r.delnolof)+"\t"+Integer.toString(r.delFS)+"\t"+Integer.toString(r.sodel));
			
			w.closew();
		}
		

	// ecdf of read coverage dbsnp sv no dbsnp
		public static void writeCovgatk(LoFWorkflowResult r, String dbfile, String nodbfile){
			
			Writer_Output w1 = new Writer_Output(dbfile);
			Writer_Output w2 = new Writer_Output(nodbfile);
			
			w1.writeFileln("varaint\tall");
			w2.writeFileln("variant\tall");
			
			Iterator<VCFLine> it = r.varlist.iterator();
			
			VCFLine l;
			while(it.hasNext()){
				
				l=it.next();
				String [] tags= l.getSample().split(":");
				String [] cov= tags[1].split(",");
				int ref=Integer.parseInt(cov[0]);
				int var=Integer.parseInt(cov[1]);
				int all=ref+var;
				
				if(l.getTool()==VCFFile.GATK){
					if(l.getID().equals(".")){
						w2.writeFileln(var+"\t"+all);
					}
					else{
						w1.writeFileln(var+"\t"+all);
					}
				}
			}
			
			w1.closew();
			w2.closew();
		}
		
	// box plot compare gatk only and gatk+pindel quality scores
		public static void writeQualgatkpindel(LoFWorkflowResult r, String qualgatk, String qualgatkpindel){
			
			
			Writer_Output w1 = new Writer_Output(qualgatk);
			Writer_Output w2 = new Writer_Output(qualgatkpindel);
			
			w1.writeFileln("qualities gatk");
			w2.writeFileln("qualities gatk+pindel");
			
			Iterator<VCFLine> it = r.varlist.iterator();
			
			VCFLine l;
			while(it.hasNext()){
				
				l=it.next();
				
				if(l.getSnpindel()!=VCFLine.SNP && l.getSnpindel()!=0){
					if(l.getTool()==VCFFile.GATK){
						w1.writeFileln(l.getQual());
					}
					else if(l.getTool()==VCFFile.GATK+VCFFile.PINDEL){
						w2.writeFileln(l.getQual());
					}
				}
			}
			
			w1.closew();
			w2.closew();
			
		}
		
	// gatk and pindel indel counts
		public static void writeIndelCounts(LoFWorkflowResult r, String counttable){

			int pindelall=0;
			int gatkall=0;
			int bothall=0;
			
			Iterator<VCFLine> it = r.varlist.iterator();
			while(it.hasNext()){
				VCFLine l = it.next();
				
				if(l.getSnpindel()== VCFLine.INS ||l.getSnpindel()== VCFLine.DEL || l.getSnpindel()== VCFLine.INDEL ){
					
					if(l.getTool()==VCFFile.GATK){
						gatkall++;
					}
					else if(l.getTool()==VCFFile.PINDEL){
						pindelall++;
					}
					else{
						bothall++;
					}
				}
			}
			
			int pindellof=0;
			int gatklof=0;
			int bothlof=0;
			
			Iterator<VCFLine> j = r.loflist.iterator();
			while(j.hasNext()){
				VCFLine l = j.next();
				
				if(l.getSnpindel()== VCFLine.INS ||l.getSnpindel()== VCFLine.DEL || l.getSnpindel()== VCFLine.INDEL ){
					
					if(l.getTool()==VCFFile.GATK){
						gatklof++;
					}
					else if(l.getTool()==VCFFile.PINDEL){
						pindellof++;
					}
					else{
						bothlof++;
					}
				}
			}
			
			Writer_Output w = new Writer_Output(counttable);
			w.writeFileln("gakt onyl indels\tpindel only indels\tgatk and pindel indels"+
					"\tgatk only fsindel\tpindel only fsindel\t pindel and gatk fsindel");
			w.writeFileln(gatkall+"\t"+pindelall+"\t"+bothall+"\t"+gatklof+"\t"+pindellof+"\t"+bothlof);
			w.closew();
			
		}

}
