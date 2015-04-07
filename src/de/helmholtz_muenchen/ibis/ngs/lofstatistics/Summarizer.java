package de.helmholtz_muenchen.ibis.ngs.lofstatistics;

import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;
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
	private HashMap<String,Integer> vep_header; //contains indices for LoF, LoF_filter etc.
	private int nrOfVEPfields;
	
	//gene_id -> sample_id -> [LoF transcripts homo, LoF transcripts hetero]
	private HashMap<String, HashMap<String,Integer []>> gene_sample_statistic;
	//gene_id -> [full LoFs, part LoFs, observed LoFs homo, observed LoFs hetero]
	private HashMap<String, Integer[]> gene_statistic;
	
	
	//sample_id -> fields of SampleStat
	private HashMap<String, SampleStat> sample_statistic;
	
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
		
		lof_statistic = new ArrayList<LoFVariant>();
		
	}
	
	public String getLoFStatistic() {
		if(lof_statistic.size()==0) {
			this.extract_LOFs();
		}
		String outfile = vcf_file.replace("vcf", "lof_statistic.txt");
		try {
			this.writeLOFStatistics(outfile);
		} catch (IOException e) {
			System.err.println(outfile+ " could not be written!");
			e.printStackTrace();
		}
		return outfile;
	}
	
	public String getGeneStatistic() {
		if(lof_statistic.size()==0) {
			this.extract_LOFs();
		}
		String outfile = vcf_file.replace("vcf","gene_statistic.txt");
		try {
			this.generateGeneStatistic();
			this.writeGeneStatistic(outfile);
		} catch (IOException e) {
			System.err.println(outfile+ " could not be written!");
			e.printStackTrace();
		}
		return outfile;
	}
	
	public String getSampleStatistic() {
		if(lof_statistic.size()==0) {
			this.extract_LOFs();
			this.generateGeneStatistic();
		}
		this.generateSampleStatistic();
		String outfile = vcf_file.replace("vcf", "sample_statistic.txt");
		try {
			this.writeSampleStatistic(outfile);
		} catch (IOException e) {
			System.err.println(outfile+" could not be written!");
		}
		return outfile;
	}
	
	public String getTranscriptStatistic() {
		if(lof_statistic.size()==0) {
			this.extract_LOFs();
		}
		this.generateTranscriptStatistic();
		String outfile = vcf_file.replace("vcf", "transcript_statistic.txt");
		try {
			this.writeTranscriptStatistic(outfile);
		} catch (IOException e) {
			System.err.println(outfile+" could not be written!");
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
				i = i.replaceFirst("CSQ=","");
				transcript_annotations = i.split(",");
				break;
			}
		}
		
		if(transcript_annotations == null) return; //nothing was annotated
		
		//fill HashMap genes
		for(String t: transcript_annotations) {
			
			annotation_fields = t.split("\\|");
			
			//work around for annotations ending with ...||
			String [] tmp = new String [nrOfVEPfields];
			for(int i = 0; i< annotation_fields.length;i++) {
				tmp[i] = annotation_fields[i];
			}
			for(int j = annotation_fields.length; j < tmp.length; j++) {
				tmp[j] = "";
			}
			annotation_fields = tmp;
			
			int confidence_index = vep_header.get("confidence");
//			System.out.println(t);
			if(annotation_fields[confidence_index].equals("")) continue;
			
			
			
			//alt_allele - ref_allele = annotation_fields[0];
			int allele_index = vep_header.get("allele");
			if(alt.contains(annotation_fields[allele_index]) || annotation_fields[allele_index].equals("-")) { //found a LoF annotation
				consequence = annotation_fields[vep_header.get("consequence")];
				gene_symbol = annotation_fields[vep_header.get("gene_symbol")];
				gene_id = annotation_fields[vep_header.get("gene_id")];
				transcript_id = annotation_fields[vep_header.get("transcript_id")];
				confidence = annotation_fields[confidence_index];
				filter = annotation_fields[vep_header.get("lof_filter")];
				lof_info = annotation_fields[vep_header.get("lof_info")];
				lof_flags = annotation_fields[vep_header.get("lof_flags")];
				
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
		
		if(genes.size()==0) return; //no LoF genes have been found

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
	
	private void extract_LOFs() {
		
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
			        	if(line.startsWith("##INFO=<ID=CSQ")) {
			        		vep_header = new HashMap<String,Integer>();
			        		String [] tmp = line.split("\\|");
			        		nrOfVEPfields = tmp.length;
			        		String part; 
			        		for(int i = 0; i<tmp.length; i++) {
			        			part = tmp[i];
			        			if(part.contains("Allele")) {
			        				vep_header.put("allele", i);
			        			} else if (part.contains("Consequence")) {
			        				vep_header.put("consequence", i);
			        			} else if (part.equals("SYMBOL")) {
			        				vep_header.put("gene_symbol",i);
			        			} else if (part.contains("Gene")) {
			        				vep_header.put("gene_id", i);
			        			} else if (part.contains("Feature")) {
			        				vep_header.put("transcript_id",i);
			        			}  else if (part.contains("LoF_filter")) {
			        				vep_header.put("lof_filter",i);
			        			} else if (part.contains("LoF_flags")) {
			        				vep_header.put("lof_flags",i);
			        			} else if (part.contains("LoF_info")) {
			        				vep_header.put("lof_info",i);
			        			} else if (part.contains("LoF")) {
			        				vep_header.put("confidence",i);
			        			}
			        		}
			        	} else if(line.startsWith("#CHR")) {
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
//			        if(counter > 50000) {
//			        	return;
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
	
	private void generateGeneStatistic() {
		//gene_id -> sample_id -> [LoF transcripts homo, LoF transcripts hetereo]
		gene_sample_statistic = new HashMap<>();
		
		HashMap<String, HashMap<String, HashSet<String>>> transcripts_hom = new HashMap<>();
		HashMap<String, HashMap<String, HashSet<String>>> transcripts_het = new HashMap<>();
		
		//gene_id -> [full LoFs, part LoFs, observed LoFs homo, observed LoFs hetero]
		gene_statistic = new HashMap<>();
		
		
		for(LoFVariant lof: lof_statistic)  {
			HashMap<String, LoFGene> gene2lof = lof.getGene_id2lofgene_fields();
			
			for(String gene_id: gene2lof.keySet()) {
				
				Integer [] stats = gene_statistic.get(gene_id);
				if(stats==null) {
					stats = new Integer[4];
					for(int i=0; i<4;i++) {
						stats[i]=0;
					}
				}
				
				HashSet<String> transcripts = gene2lof.get(gene_id).getTranscripts();
				
				//is variant full LOF
				if(transcripts.size() == gene_id2transcript_ids.get(gene_id).size()) {
					stats[0]++;
				} else {
					stats[1]++;
				}
				
				for(int i=0; i<lof.getGenotypes().size(); i++) {
					String gt = lof.getGenotypes().get(i);
					
					//homozygote affected transcripts
					if(gt.charAt(0)==gt.charAt(2) && gt.charAt(0)=='1') {
						stats[2]++; //LoF is homozygous
						HashMap<String, HashSet<String>> sample2transcripts = transcripts_hom.get(gene_id);
						if(sample2transcripts==null) {
							sample2transcripts = new HashMap<>();
						}
						HashSet<String> tmp_set = sample2transcripts.get(sample_ids.get(i));
						if(tmp_set == null) {
							tmp_set = new HashSet<>();
						}
						
						tmp_set.addAll(transcripts);
						sample2transcripts.put(sample_ids.get(i), tmp_set);
						transcripts_hom.put(gene_id, sample2transcripts);
					} else if (gt.charAt(0)=='1' || gt.charAt(2)=='1') {
						stats[3]++; //LoF is heterozygous
						HashMap<String, HashSet<String>> sample2transcripts = transcripts_het.get(gene_id);
						if(sample2transcripts==null) {
							sample2transcripts = new HashMap<>();
						}
						HashSet<String> tmp_set = sample2transcripts.get(sample_ids.get(i));
						if(tmp_set == null) {
							tmp_set = new HashSet<>();
						}
						
						tmp_set.addAll(transcripts);
						sample2transcripts.put(sample_ids.get(i), tmp_set);
						transcripts_het.put(gene_id, sample2transcripts);
					}
				}
				gene_statistic.put(gene_id, stats);
			}
		}
		
		int hom_lof_nr, het_lof_nr;
		for(String gene: gene_id2gene_symbol.keySet()) {
			for(String sample: sample_ids) {
				HashSet<String> hom_lofs = null;
				HashSet<String> het_lofs = null;
				
				if(transcripts_hom.get(gene)!=null) {
					hom_lofs = transcripts_hom.get(gene).get(sample);
				}
				if(transcripts_het.get(gene)!=null) {
					het_lofs = transcripts_het.get(gene).get(sample);
				}
				
				if(hom_lofs==null) hom_lofs = new HashSet<>();
				if(het_lofs==null) het_lofs = new HashSet<>();
				
				het_lofs.removeAll(hom_lofs);
				hom_lof_nr = hom_lofs.size();
				het_lof_nr = het_lofs.size();
				HashMap<String, Integer []> sample2nr = gene_sample_statistic.get(gene);
				if(sample2nr == null) sample2nr = new HashMap<String, Integer[]>();
				sample2nr.put(sample, new Integer[]{hom_lof_nr,het_lof_nr});
				gene_sample_statistic.put(gene, sample2nr);
			}
		}
	}
	
	private void generateSampleStatistic() {
		this.sample_statistic = new HashMap<>();
		ArrayList<String> completeLOFgene, partLOFgene;
		HashMap<String, Integer> fullLOFs = new HashMap<>();
		HashMap<String, Integer> partLOFs = new HashMap<>();
		boolean isFull = false;
		
		//fill fullLOFs and partLOFs
		for(LoFVariant lof:lof_statistic) {
			HashMap<String, LoFGene> lof_genes = lof.getGene_id2lofgene_fields();
			isFull = false;
			for(Entry<String, LoFGene> gene:lof_genes.entrySet()) {
				int all_trans = gene_id2transcript_ids.get(gene.getKey()).size();
				int lof_trans = gene.getValue().getNrOfTranscripts();
				if(all_trans == lof_trans) {
					isFull = true;
				}
			}
			for(int i = 0; i< lof.getGenotypes().size(); i++) {
				String gt = lof.getGenotypes().get(i);
				String sample = sample_ids.get(i);
				if(gt.charAt(0)=='1' || gt.charAt(2)=='1') {
					if(isFull) {
						Integer t = fullLOFs.get(sample);
						if(t==null) {
							t = 0;
						}
						t++;
						fullLOFs.put(sample, t);
					} else {
						Integer x = partLOFs.get(sample);
						if(x == null) {
							x = 0;
						}
						x++;
						partLOFs.put(sample, x);
					}
				}
			}
		}
		
		//get complete LOF genes
		for(String sample: sample_ids) {
			
			completeLOFgene = new ArrayList<>();
			partLOFgene = new ArrayList<>();
			
			for(String gene: gene_id2gene_symbol.keySet()) {
				int all_trans = gene_id2transcript_ids.get(gene).size();
				int lof_hom = gene_sample_statistic.get(gene).get(sample)[0];
				if(all_trans == lof_hom) {
					completeLOFgene.add(gene);
				} else {
					partLOFgene.add(gene);
				}
			}
			sample_statistic.put(sample, new SampleStat(fullLOFs.get(sample),partLOFs.get(sample), completeLOFgene, partLOFgene));
		}
	}
	
	private void generateTranscriptStatistic() {
		transcript_statistic = new HashMap<>();
		for(LoFVariant lofvar: lof_statistic) {
			HashSet<String> transcripts = new HashSet<>();
			HashMap<String, LoFGene> lofgene = lofvar.getGene_id2lofgene_fields();
			int obs_hom = 0;
			int obs_het = 0;
			for(String g: lofgene.keySet()) {
				transcripts.addAll(lofgene.get(g).getTranscripts());
			}
			for(String gt: lofvar.getGenotypes()) {
				if(gt.charAt(0)=='1' && gt.charAt(2)==gt.charAt(0)) {
					obs_hom ++;
				} else if(gt.charAt(0)=='1' || gt.charAt(2)=='1') {
					obs_het ++;
				}
			}
			for(String t: transcripts) {
				Integer[] stats = transcript_statistic.get(t);
				if(stats==null) {
					stats = new Integer[3];
					for (int i = 0; i< 3; i++) {
						stats[i]= 0;
					}
				}
				stats[0]++;
				stats[1]+= obs_hom;
				stats[2]+= obs_het;
				transcript_statistic.put(t, stats);
			}
		}
	}
	
	private void writeLOFStatistics(String outfile) throws IOException{
		
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
	
	private void writeGeneStatistic(String outfile) throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(outfile));
		
		bw.write("gene_id\tgene_symbol\tfull\tpartial\tobs_hom\tobs_het\tall_trans");
		for(String s: sample_ids) {
			bw.write("\t"+s);
		}
		bw.newLine();
		for(String gene: gene_id2gene_symbol.keySet()) {
			bw.write(gene+"\t"+gene_id2gene_symbol.get(gene));
			for(int i: gene_statistic.get(gene)){
				bw.write("\t"+i);
			}
			bw.write("\t"+gene_id2transcript_ids.get(gene).size());
			for(String sample: sample_ids) {
				int hom_lof = gene_sample_statistic.get(gene).get(sample)[0];
				int het_lof = gene_sample_statistic.get(gene).get(sample)[1];
				bw.write("\t"+hom_lof+","+het_lof);
			}
			bw.newLine();
		}
		
		bw.close();
	}
	
	private void writeSampleStatistic(String outfile) throws IOException {
		
		BufferedWriter bw = new BufferedWriter(new FileWriter(outfile));
		bw.write("sample_id\tfull\tpartial\tknocked_out\tcompleteLOFgenes");
		bw.newLine();
		for(String s: sample_statistic.keySet()) {
			SampleStat stat = sample_statistic.get(s);
			ArrayList<String> completes = stat.getComplete_LOF_genes();
			bw.write(s+"\t"+stat.getFullLOFs()+"\t"+stat.getPartLOFs()+"\t"+completes.size()+"\t");
			if(completes.size() >= 1) {
				String gene = completes.remove(0);
				bw.write(gene +":"+gene_id2gene_symbol.get(gene));
			}
			
			for(String u: completes) {
				bw.write(","+u+":"+gene_id2gene_symbol.get(u));
			}
			bw.newLine();
		}
		bw.close();
	}
	
	private void writeTranscriptStatistic(String outfile) throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(outfile));
		
		bw.write("transcript_id\tgene_id\tgene_symbol\tlofs\tobs_hom\tobs_het");
		bw.newLine();
		
		for(String t: transcript_statistic.keySet()) {
			String gene_id = transcript_id2gene_id.get(t);
			String gene_symbol = gene_id2gene_symbol.get(gene_id);
			Integer [] stats = transcript_statistic.get(t);
			bw.write(t+"\t"+gene_id+"\t"+gene_symbol+"\t"+stats[0]+"\t"+stats[1]+"\t"+stats[2]);
			bw.newLine();
		}
		
		bw.close();
	}
	
	public static void main (String [] args) {
		String vcf_file;
		String cds_file;
		
		/**@ helmholtz**/
		vcf_file = "/home/ibis/tim.jeske/KORAdata/analysis_ready/20150401_loftee_q/variant_effect_LoF.vcf";
		cds_file = "/home/ibis/tim.jeske/LOFTEE/Homo_sapiens.GRCh37.75.cds.all.fa";
		
		Summarizer my = new Summarizer("LOFTEE",vcf_file,cds_file);
//		my.getLoFStatistic();
//		my.getGeneStatistic();
//		my.getSampleStatistic();
		my.getTranscriptStatistic();
	}
	
}