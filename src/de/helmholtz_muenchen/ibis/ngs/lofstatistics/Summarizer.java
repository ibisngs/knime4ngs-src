package de.helmholtz_muenchen.ibis.ngs.lofstatistics;

import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Scanner;

public class Summarizer {

	private String annotation;
	private String vcf_file;
	private String cds_file;
	
	//internal data containers
	private HashMap<String,String> gene_id2gene_symbol; //filled while calculation!
	private HashMap<String,String> transcript_id2gene_id;
	private HashMap<String,ArrayList<String>> gene_id2transcript_ids;
	private ArrayList<String> sample_ids;
	
	//statistics
	
	//main statistic
	private ArrayList<LoFVariant> lof_statistic;
	private ArrayList<String> additional_titles;
	private ArrayList<ArrayList<HashSet<String>>> additional_fields;
	
	//gene_id -> sample_id -> [LoF transcripts homo, LoF transcripts hetereo]
	private HashMap<String, HashMap<String,Integer []>> gene_sample_statistic;
	
	//sample_id -> fields of SampleStat
	private HashMap<String, SampleStat> sample_statistic;
	
	//gene_id -> [full LoFs, part LoFs, observed LoFs homo, observed LoFs hetero]
	private HashMap<String, Integer[]> gene_statistic;
	
	//transcript_id -> [LoFs, observed LoFs homo, observed LoFs hetero]
	private HashMap<String, Integer[]> transcript_statistic;
	
	/**
	 * 
	 * @param annotation LOFTEE, VAT or others
	 * @param vcf_file file containing annotation
	 * @param cds_file file containing gene id for each transcript id, can be found
	 *                 here: http://www.ensembl.org/info/data/ftp/index.html
	 */
	public Summarizer(String annotation, String vcf_file, String cds_file) {
		this.annotation = annotation;
		this.vcf_file = vcf_file;
		this.cds_file = cds_file;
		try {
			System.out.println("Reading CDS file...");
			this.readCDSFile();
			System.out.println("CDS file read.");
		} catch (IOException e) {
			System.err.println("Reading CDS file failed. LoF effect (full or partial) cannot be calculated.");
			e.printStackTrace();
		}
		
		additional_titles = new ArrayList<String>();
		
		//TODO add another annotation
		if(annotation.equals("LOFTEE")) {
			additional_titles.add("confidence");
			additional_titles.add("failed_filters");
			additional_titles.add("lof_flags");
			additional_titles.add("lof_info");
		}
		
		additional_fields = new ArrayList<ArrayList<HashSet<String>>>();
		
		lof_statistic = new ArrayList<LoFVariant>();
		
	}
	
	public String getLoFStatistic() {
		this.extract_LOFs();
		String outfile = vcf_file.replace("vcf", "lofstatistic.vcf");
		try {
			this.writeLOFStatistics(outfile);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return outfile;
	}
	
	/**
	 * reads CDS file and fills transcript_id2gene_id and gene_id2transcript_ids
	 * @throws IOException
	 */
	private void readCDSFile() throws IOException {
		
		transcript_id2gene_id = new HashMap<String,String>();
		gene_id2transcript_ids = new HashMap<String, ArrayList<String>>();
		
		String [] fields;
		String transcript_id, gene_id;
		
		FileInputStream inputStream = null;
		Scanner sc = null;
		try {
		    inputStream = new FileInputStream(this.cds_file);
		    sc = new Scanner(inputStream, "UTF-8");
		    while (sc.hasNextLine()) {
		        String line = sc.nextLine();
		        if(line.startsWith(">")) {
		        	fields = line.split("\\s");
		        	transcript_id = fields[0].replaceFirst(">","");
		        	gene_id = fields[3].split(":")[1];
		        	transcript_id2gene_id.put(transcript_id, gene_id);
		        	
		        	if(gene_id2transcript_ids.containsKey(gene_id)) {
		        		gene_id2transcript_ids.get(gene_id).add(transcript_id);
		        	} else {
		        		ArrayList<String> transcript_ids = new ArrayList<String>();
		        		transcript_ids.add(transcript_id);
		        		gene_id2transcript_ids.put(gene_id, transcript_ids);
		        	}
		        }
		    }
		    // note that Scanner suppresses exceptions
		    if (sc.ioException() != null) {
		        throw sc.ioException();
		    }
		} finally {
		    if (inputStream != null) {
		        inputStream.close();
		    }
		    if (sc != null) {
		        sc.close();
		    }
		}
	}
	
	private String replaceVATbySO(String type) {
		String result = null;
		if(type.equals("spliceOverlap")) {
			result = "splice_site_variant";
		} else if(type.equals("prematureStop")){
			result = "stop_gained";
		} else if(type.equals("insertionFS") || type.equals("deletionFS")) { 
			result = "frameshift_variant";
		}
		return result;
	}
	
	private void getLoFVariantfromVAT(String chr, String pos, String ref, String id, String alt, int GT_index,String[] infos, ArrayList<String> genotypes) {
		
		HashMap<String,LoFGene> genes = new HashMap<String,LoFGene>();
		
		int nrOfTitles = additional_titles.size();
		
		//init additional fields
		
		int observed_homo = 0;
		int observed_hetero = 0;
		String gene_id= "";
		String consequence="";
		String gene_symbol;
		
		String [] allele_annotations = null;
		String [] annotation_fields = null;
		
		for(String i: infos) {
			if(i.startsWith("VA")) {
				i = i.replaceFirst("VA=","");
				allele_annotations = i.split(",");
				break;
			}
		}
		
		if(allele_annotations == null) return; //nothing was annotated
		
		for(String a: allele_annotations) {
			annotation_fields = a.split(":");
			if(annotation_fields[0].equals(GT_index+"")) {
				consequence = replaceVATbySO(annotation_fields[4]);
				if(consequence==null) return; //no LoF variant annotated
				gene_id = annotation_fields[2];
				if(gene_id.contains(".")) {
					gene_id = gene_id.split("\\.")[0];
				}
				gene_symbol = annotation_fields[1];
				gene_id2gene_symbol.put(gene_id, gene_symbol);
				
				LoFGene g;
				if(genes.containsKey(gene_id)) {
					g = genes.get(gene_id);
				} else {
					g = new LoFGene(nrOfTitles);
				}
				g.addConsequence(consequence);
				
				for(String s : annotation_fields) {
					if(s.startsWith("ENST")) {
						if(s.contains(".")) {
							g.addTranscript(s.split("\\.")[0]);
						} else {
							g.addTranscript(s);
						}
					}
				}
				genes.put(gene_id, g);
			}
		}
		
		ArrayList<String> adapted_genotypes = getAdaptedGenotypes(genotypes, GT_index);
		for(String str: adapted_genotypes) {
			if(str.charAt(0)==str.charAt(2) && str.charAt(0)=='1') {
				observed_homo++;
			} else if (str.charAt(0)=='1' || str.charAt(2)=='1') {
				observed_hetero++;
			}
		}
		lof_statistic.add(new LoFVariant(chr, pos, ref, alt, id, observed_homo, observed_hetero, genes, adapted_genotypes));
	}
	
	/**
	 * fills lof_statistic and gene_id2gene_symbol
	 */
	private void getLoFVariantfromLOFTEE(String chr, String pos, String ref, String id, String alt, int GT_index,String[] infos, ArrayList<String> genotypes) {
		
		
		HashMap<String,LoFGene> genes = new HashMap<String,LoFGene>();
		
		int nrOfTitles = additional_titles.size();
		
		//init additional fields
		
		int observed_homo = 0;
		int observed_hetero = 0;
		String [] transcript_annotations=null;
		String [] annotation_fields;
		String gene_id= "";
		String consequence, transcript_id, confidence, filter;
		String gene_symbol;
		String lof_info = "";
		String lof_flags = "";
		
		for(String i:infos) {
			if(i.startsWith("CSQ")) {
				if(!i.contains("|HC|") && !i.contains("|LC|")) return; //no LoF variants were annotated
				i = i.replaceFirst("CSQ=","");
				transcript_annotations = i.split(",");
				break;
			}
		}
		
		if(transcript_annotations == null) return; //nothing was annotated
		
		//fill HashMap genes
		for(String t: transcript_annotations) {
			
			if(!t.contains("|HC|") && !t.contains("|LC|")) continue;
//			System.out.println(t);
			annotation_fields = t.split("\\|");
			
			//alt_allele - ref_allele = annotation_fields[0];
			if(alt.contains(annotation_fields[0]) || annotation_fields[0].equals("-")) { //found a LoF annotation
				consequence = annotation_fields[1];
				gene_symbol = annotation_fields[3];
				gene_id = annotation_fields[4];
				transcript_id = annotation_fields[6];
				confidence = annotation_fields[22];
				filter = annotation_fields[23];
				
				if(annotation_fields.length>24) {
					lof_info = annotation_fields[24];
				}
				if(annotation_fields.length>25) {
					lof_flags = annotation_fields[25];
				}
				LoFGene g;
				if(genes.containsKey(gene_id)) {
					g = genes.get(gene_id);
				} else {
					g = new LoFGene(nrOfTitles);
				}
				g.addConsequence(consequence);
				g.addTranscript(transcript_id);
				
				g.addInfo(confidence, additional_titles.indexOf("confidence"));
				g.addInfo(filter, additional_titles.indexOf("failed_filters"));
				g.addInfo(lof_info, additional_titles.indexOf("lof_info"));
				g.addInfo(lof_flags, additional_titles.indexOf("lof_flags"));
				
				genes.put(gene_id, g);
//				System.out.println(genes.size());
				gene_id2gene_symbol.put(gene_id, gene_symbol);
			}
		}

		ArrayList<String> adapted_genotypes = getAdaptedGenotypes(genotypes, GT_index);
		for(String str: adapted_genotypes) {
			if(str.charAt(0)==str.charAt(2) && str.charAt(0)=='1') {
				observed_homo++;
			} else if (str.charAt(0)=='1' || str.charAt(2)=='1') {
				observed_hetero++;
			}
		}
		lof_statistic.add(new LoFVariant(chr, pos, ref, alt, id, observed_homo, observed_hetero, genes, adapted_genotypes));
	}
	
	
	private ArrayList<String> getAdaptedGenotypes(ArrayList<String> genotypes, int GT_index) {
		ArrayList<String> adapted_genotypes = new ArrayList<String>();
		for(String gt:genotypes) {
			String [] genotype;
			String separator;
			if(gt.contains("/")) {
				separator = "/";
				genotype = gt.split("/");
			} else {
				separator = "|";
				genotype = gt.split("\\|");
			}
			String first = genotype[0];
			String second = genotype[1];
			String first_a = ".";
			String second_a = ".";
			if(first.equals(second) && first.equals(GT_index+"")) {
				first_a = "1";
				second_a = "1";
			} else if (first.equals(GT_index+"")) {
				first_a = "1";
				if(second.equals("0")) {
					second_a = "0";
				}
			} else if (second.equals(GT_index+"")) {
				second_a = "1";
				if(first.equals("0")) {
					first_a = "0";
				}
			} else {
				first_a = first;
				second_a = second;
			}
			String adapted = first_a+separator+second_a;
//			if(!adapted.equals(gt)) {
//				System.out.println("Changed genotype: "+gt+" to "+adapted);
//			}
			adapted_genotypes.add(adapted);
		}
		return adapted_genotypes;
	}
	
	private  ArrayList<String> getGenotypes(String [] fields) {
		ArrayList<String> result = new ArrayList<String>();
		String gt = "";
		
		for(int i=9;i<fields.length;i++) {
			gt = fields[i].split(":")[0];
			result.add(gt);
		}
		return result;
	}
	
	public void extract_LOFs() {
		
		gene_id2gene_symbol = new HashMap<String,String>();
		
		try {
			System.out.println("Processing VCF file...");
			int counter = 0;
			String [] fields;
			String chr, pos, ref;
			String [] ids, alt_alleles, infos;
			ArrayList<String> genotypes;
			FileInputStream inputStream = null;
			Scanner sc = null;
			try {
			    inputStream = new FileInputStream(vcf_file);
			    sc = new Scanner(inputStream, "UTF-8");
			    while (sc.hasNextLine()) {
			    	if(counter%10000==0) {
			    		System.out.println("Processed "+counter+" lines...");
			    	}
			        String line = sc.nextLine();
			        if(line.startsWith("#"))  {
			        	if(line.startsWith("#CHR")) {
			        		String [] header = line.split("\t");
			        		sample_ids = new ArrayList<String>();
			        		for(int i = 9; i< header.length; i++) {
			        			sample_ids.add(header[i]);
			        		}
			        	}
			        } else {
			        	fields = line.split("\t");
			        	chr = fields[0];
			        	pos = fields[1];
			        	ref = fields[3];
			        	
			        	ids = fields[2].split(";");
			        	alt_alleles = fields[4].split(",");
			        	infos = fields[7].split(";");
			        	genotypes = getGenotypes(fields);
			        	
			        	for(int i = 0; i<alt_alleles.length;i++) {
			        		String id = ".";
			        		if(ids.length>i) {
			        			id = ids[i];
			        		}
//			        		TODO implement some method to retrieve id from dbSNP
//			        		if(id.equals(".")) {
//			        			getId(chr, pos, ref, alt_alleles[i]);
//			        		}
			        		if(annotation.equals("LOFTEE")) {
			        			getLoFVariantfromLOFTEE(chr, pos, ref, id, alt_alleles[i], i+1,infos, genotypes);
			        		} else if (annotation.equals("VAT")) {
			        			getLoFVariantfromVAT(chr, pos, ref, id, alt_alleles[i], i+1,infos, genotypes);
			        		}
			        	}
			        }
			        counter++;
			        /**only for testing**/
//			        if(counter > 20000) {
//			        	this.writeLOFStatistics();
//			        	System.exit(0);
//			        }
			    }
			    System.out.println("VCF file processed.");
			    // note that Scanner suppresses exceptions
			    if (sc.ioException() != null) {
			        throw sc.ioException();
			    }
			} finally {
			    if (inputStream != null) {
			        inputStream.close();
			    }
			    if (sc != null) {
			        sc.close();
			    }
			}
			
		} catch (IOException e) {
			System.err.println("VCF file couldn't be read!");
			e.printStackTrace();
		}
	}
	
	public  void writeLOFStatistics(String outfile) throws IOException{
		
		BufferedWriter bw = new BufferedWriter(new FileWriter(outfile));
		String [] header = {"chr", "pos", "rsId", "ref_allele", "alt_allele", "gene_id", "gene_symbol", "effect", "consequence", "lof_trans", "all_trans", "obs_hom", "obs_het"};
		
		//write header
		for(String s:header) {
			bw.write(s+"\t");
		}
		for(String s: additional_titles) {
			bw.write(s+"\t");
		}
		bw.write(sample_ids.remove(0));
		for(String s:sample_ids) {
			bw.write("\t"+s);
		}
		bw.newLine();
		
		//write lines
		for(LoFVariant l:lof_statistic) {
			HashMap<String, LoFGene> map = l.getGene_id2lofgene_fields();
			LoFGene lofgene;
			String gene_ids = "";
			String genes = "";
			String effects = "";
			String lof_trans = "";
			String all_trans = "";
			String consequences = "";
			String [] add_infos = new String[additional_titles.size()];
			for(int i = 0; i< additional_titles.size(); i++) {
				add_infos[i] = "";
			}
			
			for(String s:map.keySet()) {
				lofgene = map.get(s);
				gene_ids += ","+s;
				genes += ","+ gene_id2gene_symbol.get(s);
				String nrOfTrans = "-";
				if(gene_id2transcript_ids.containsKey(s)) {
					nrOfTrans = gene_id2transcript_ids.get(s).size()+"";	
				} else {
					System.err.println("According to the CDS file gene "+s+" has no protein coding transcripts.");
				}
				all_trans += "," + nrOfTrans;
				lof_trans += "," + lofgene.getNrOfTranscripts();
				if(nrOfTrans.equals(lofgene.getNrOfTranscripts()+"")) {
					effects += ",full";
				} else if (nrOfTrans.equals("-")){
					effects += ",unknown";
				} else {
					effects += ",partial";
				}
				
				for(String c: lofgene.getConsequences()) {
					consequences += "," + c;
				}

				for(int i=0; i < additional_titles.size(); i++) {
					HashSet<String> infos = lofgene.getInfo(i);
					for(String str:infos) {
						if(!str.equals("")) {
							add_infos[i]+=","+str;
						}
					}
				}
			}
			gene_ids = gene_ids.replaceFirst(",","");
			genes = genes.replaceFirst(",","");
			effects = effects.replaceFirst(",", "");
			consequences = consequences.replaceFirst(",", "");
			lof_trans = lof_trans.replaceFirst(",","");
			all_trans = all_trans.replaceFirst(",","");
			
			bw.write(l.getChr()+"\t"+l.getPos()+"\t"+l.getRsId()+"\t"+l.getRef_allele()+"\t"+l.getAlt_allele()+"\t"+gene_ids+"\t"+genes+"\t"+effects+"\t"+consequences+"\t"+lof_trans+"\t"+all_trans+"\t"+l.getObserved_homo()+"\t"+l.getObserved_hetero()+"\t");
			
			//write additional fields
			for(String s: add_infos) {
				s = s.replaceFirst(",","");
				bw.write(s+"\t");
			}
			
			ArrayList<String> genotypes = l.getGenotypes();
			bw.write(genotypes.remove(0));
			for(String gt: genotypes) {
				bw.write("\t"+gt);
			}
			bw.newLine();
		}
		
		bw.close();
	}
	
	public static void main (String [] args) {
		String vcf_file;
		String cds_file;
		
		/**@ home**/
//		vcf_file = "/home/tim/LOF_Project_local/files/1000genomes/chr22/variant_effect_LoF.vcf";
//		cds_file = "/home/tim/LOF_Project_local/files/Homo_sapiens.GRCh38.cds.all.fa";
		
		/**@ helmholtz**/
//		vcf_file = "/home/ibis/tim.jeske/KORAdata/analysis_ready/20150325_loftee_q/variant_effect_LoF.vcf";
		vcf_file = "/home/ibis/tim.jeske/KORAdata/analysis_ready/20150319_vat_k/analysis_ready.diabetes.filtered.haplotypecaller.INDEL.vat.vcf";
		cds_file = "/home/ibis/tim.jeske/LOFTEE/Homo_sapiens.GRCh38.cds.all.fa";
		
		
		Summarizer my = new Summarizer("VAT",vcf_file,cds_file);
		my.getLoFStatistic();
	}
	
}