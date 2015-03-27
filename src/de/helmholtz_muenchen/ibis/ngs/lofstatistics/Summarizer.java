package de.helmholtz_muenchen.ibis.ngs.lofstatistics;

import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Scanner;

public class Summarizer {

	//internal data containers
	private static HashMap<String,String> gene_id2gene_symbol; //filled while calculation!
	private static HashMap<String,String> transcript_id2gene_id;
	private static HashMap<String,ArrayList<String>> gene_id2transcript_ids;
	private static ArrayList<String> sample_ids;
	
	//statistics
	private static ArrayList<LoFVariant> lof_statistic;
	
	//gene_id -> sample_id -> [LoF transcripts homo, LoF transcripts hetereo]
	private static HashMap<String, HashMap<String,Integer []>> gene_sample_statistic;
	
	//sample_id -> fields of SampleStat
	private static HashMap<String, SampleStat> sample_statistic;
	
	//gene_id -> [full LoFs, part LoFs, observed LoFs homo, observed LoFs hetero]
	private static HashMap<String, Integer[]> gene_statistic;
	
	//transcript_id -> [LoFs, observed LoFs homo, observed LoFs hetero]
	private static HashMap<String, Integer[]> transcript_statistic;
	
	private static void init() {
		gene_id2gene_symbol = new HashMap<String,String>();
		transcript_id2gene_id = new HashMap<String,String>();
		gene_id2transcript_ids = new HashMap<String,ArrayList<String>>();
		sample_ids = new ArrayList<String>();
		lof_statistic = new ArrayList<LoFVariant>();
		gene_sample_statistic = new HashMap<String, HashMap<String,Integer[]>>();
		sample_statistic = new HashMap<String,SampleStat>();
		gene_statistic = new HashMap<String, Integer[]>();
		transcript_statistic = new HashMap<String, Integer[]>();
	}
	
	private static void readCDSFile(String cds_file) throws IOException {
		
		String [] fields;
		String transcript_id, gene_id;
		
		FileInputStream inputStream = null;
		Scanner sc = null;
		try {
		    inputStream = new FileInputStream(cds_file);
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
	
	private static LoFVariant getLoFVariantfromLOFTEE(String chr, String pos, String ref, String id, String alt, int GT_index,String[] infos, ArrayList<String> genotypes) {
		HashMap<String,LoFGene> genes = new HashMap<String,LoFGene>();
		ArrayList<String> adapted_genotypes = new ArrayList<String>();
		int observed_homo = 0;
		int observed_hetero = 0;
		String [] transcript_annotations=null;
		String [] annotation_fields;
		String gene_id, consequence, transcript_id, confidence, filter, lof_info;
		String gene_symbol;
		String lof_flags ="";
		
		for(String i:infos) {
			if(i.startsWith("CSQ")) {
				if(!i.contains("|HC|") && !i.contains("|LC|")) return null;
				i = i.replaceFirst("CSQ=","");
				transcript_annotations = i.split(",");
				break;
			}
		}
		
		if(transcript_annotations == null) return null; //nothing was annotated
		
		//fill HashMap genes
		for(String t: transcript_annotations) {
			
			if(!t.contains("|HC|") && !t.contains("|LC|")) continue;
			System.out.println(t);
			annotation_fields = t.split("\\|");
			
			/**TODO:works always?**/
			if(alt.contains(annotation_fields[0]) || annotation_fields[0].equals("-")) { //found a LoF annotation
				consequence = annotation_fields[1];
				gene_symbol = annotation_fields[3];
				gene_id = annotation_fields[4];
				transcript_id = annotation_fields[6];
				confidence = annotation_fields[22];
				filter = annotation_fields[23];
				lof_info = annotation_fields[24];
				if(annotation_fields.length==26) {
					lof_flags = annotation_fields[25];
				}
				LoFGene g;
				if(genes.containsKey(gene_id)) {
					g = genes.get(gene_id);
				} else {
					g = new LoFGene();
				}
				g.addConsequence(consequence);
				g.addTranscript(transcript_id);
				g.addConfidences(confidence);
				g.addFilter(filter);
				g.addInfo(lof_info);
				g.addFlag(lof_flags);
				genes.put(gene_id, g);
				System.out.println(genes.size());
				gene_id2gene_symbol.put(gene_id, gene_symbol);
			}
		}
		
		//calculate observed_home and observed_hetero and adapt genotype list
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
				observed_homo++;
				first_a = "1";
				second_a = "1";
			} else if (first.equals(GT_index+"")) {
				observed_hetero++;
				first_a = "1";
				if(second.equals("0")) {
					second_a = "0";
				}
			} else if (second.equals(GT_index+"")) {
				observed_hetero++;
				second_a = "1";
				if(first.equals("0")) {
					first_a = "0";
				}
			} else {
				first_a = first;
				second_a = second;
			}
			String adapted = first_a+separator+second_a;
			if(!adapted.equals(gt)) {
				System.out.println("Changed genotype: "+gt+" to "+adapted);
			}
			adapted_genotypes.add(adapted);
		}
		
		return new LoFVariant(chr, pos, ref, alt, id, observed_homo, observed_hetero, genes, adapted_genotypes);
	}
	
	private static ArrayList<String> getGenotypes(String [] fields) {
		ArrayList<String> result = new ArrayList<String>();
		String gt = "";
		
		for(int i=9;i<fields.length;i++) {
			gt = fields[i].split(":")[0];
			result.add(gt);
		}
		
		return result;
	}
	
	public static ArrayList<LoFVariant> extract_LOFs(String vcf_file, String cds_file, String annotation) {
		init();
		
		ArrayList<LoFVariant> lofs = new ArrayList<LoFVariant>();
		
		try {
			System.out.println("Reading CDS file...");
			readCDSFile(cds_file);
			System.out.println("CDS file read.");
			
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
			    	if(counter%1000==0) {
			    		System.out.println("Processed "+counter+" lines...");
			    	}
			        String line = sc.nextLine();
			        if(line.startsWith("#"))  {
			        	if(line.startsWith("#CHR")) {
			        		String [] header = line.split("\t");
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
			        			LoFVariant var = getLoFVariantfromLOFTEE(chr, pos, ref, id, alt_alleles[i], i+1,infos, genotypes);
			        			if(var != null) {
			        				lofs.add(var);
			        			}
			        		}
			        	}
			        }
			        counter++;
			        /**only for testing**/
			        if(counter > 50000) {
			        	writeLOFStatistic(lofs,"/home/tim/LOF_Project_local/files/1000genomes/chr22/variant_effect_LoF.lofstatistic.vcf");
			        	System.exit(0);
			        }
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
			System.err.println("CDS file couldn't be read!");
			e.printStackTrace();
		}
		
		return lofs;
	}
	
	public static void writeLOFStatistic(ArrayList<LoFVariant> lofs, String outfile) throws IOException{
		BufferedWriter bw = new BufferedWriter(new FileWriter(outfile));
		String [] header = {"chr", "pos", "rsId", "ref_allele", "alt_allele", "gene_id", "gene_symbol", "effect", "consequence", "lof_trans", "all_trans", "confidence", "failed_filters", "lof_flags", "lof_info", "obs_hom", "obs_het"};
		
		//write header
		for(String s:header) {
			bw.write(s+"\t");
		}
		bw.write(sample_ids.remove(0));
		for(String s:sample_ids) {
			bw.write("\t"+s);
		}
		bw.newLine();
		
		//write lines
		for(LoFVariant l:lofs) {
			HashMap<String, LoFGene> map = l.getGene_id2lofgene_fields();
			LoFGene lofgene;
			String gene_ids = "";
			String genes = "";
			String effects = "";
			String lof_trans = "";
			String all_trans = "";
			String consequences = "";
			String confidences = "";
			String failedfilters = "";
			String lofflag = "";
			String lofinfo = "";
			
			for(String s:map.keySet()) {
				lofgene = map.get(s);
				gene_ids += ","+s;
				genes += ","+ gene_id2gene_symbol.get(s);
				String nrOfTrans = "-";
				if(gene_id2transcript_ids.containsKey(s)) {
					nrOfTrans = gene_id2transcript_ids.get(s).size()+"";	
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
				
				for(String c: lofgene.getConfidences()) {
					confidences += "," + c;
				}
				
				for(String f: lofgene.getFilter()) {
					failedfilters += "," + f;
				}
				
				for(String f: lofgene.getFlags()) {
					lofflag += "," + f;
				}
				
				for(String f: lofgene.getInfos()) {
					lofinfo += "," + f;
				}
				
			}
			gene_ids = gene_ids.replaceFirst(",","");
			genes = genes.replaceFirst(",","");
			effects = effects.replaceFirst(",", "");
			consequences = consequences.replaceFirst(",", "");
			lof_trans = lof_trans.replaceFirst(",","");
			all_trans = all_trans.replaceFirst(",","");
			confidences = confidences.replaceFirst(",", "");
			failedfilters = failedfilters.replaceFirst(",", "");
			lofflag = lofflag.replaceFirst(",","");
			lofinfo = lofinfo.replaceFirst(",", "");
			
			bw.write(l.getChr()+"\t"+l.getPos()+"\t"+l.getRsId()+"\t"+l.getRef_allele()+"\t"+l.getAlt_allele()+"\t"+gene_ids+"\t"+genes+"\t"+effects+"\t"+consequences+"\t"+lof_trans+"\t"+all_trans+"\t"+confidences+"\t"+failedfilters+"\t"+lofflag+"\t"+lofinfo+"\t"+l.getObserved_homo()+"\t"+l.getObserved_hetero()+"\t");
			
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
		String vcf_file = "/home/tim/LOF_Project_local/files/1000genomes/chr22/variant_effect_LoF.vcf";
		String cds_file = "/home/tim/LOF_Project_local/files/Homo_sapiens.GRCh38.cds.all.fa";
		ArrayList<LoFVariant> mylist = extract_LOFs(vcf_file, cds_file, "LOFTEE");
		try {
			writeLOFStatistic(mylist,vcf_file.replace("vcf", "lofstatistic.vcf"));
		} catch (IOException e) {
			System.err.println("lofstatistic could not be written out in file: "+vcf_file.replace("vcf", "lofstatistic.vcf"));
			e.printStackTrace();
		}
	}
	
}