package de.helmholtz_muenchen.ibis.ngs.lofsummary;

import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;
import java.util.Scanner;

import org.knime.core.node.NodeLogger;

public abstract class Summarizer {
	
	String vcf_file;
	String cds_file;
	
	//internal data containers
	HashMap<String,String> gene_id2gene_symbol; //filled while calculation!
	HashMap<String,String> transcript_id2gene_id;
	HashMap<String,ArrayList<String>> gene_id2transcript_ids;
	ArrayList<String> sample_ids;
	
	//summaries
	
	//variant summary
	ArrayList<LoFVariant> lof_statistic;
	ArrayList<String> additional_titles;
	HashMap<String,Integer> vep_header; //contains indices for LoF, LoF_filter etc.
	int nrOfVEPfields;
	
	//gene_id -> sample_id -> [LoF transcripts homo, LoF transcripts hetero]
	HashMap<String, HashMap<String,Integer []>> gene_sample_statistic;
	
	//gene_id -> sample_id -> is_unaffected
	HashMap<String, HashMap<String,Boolean>> gene_sample_unaffected;
	//gene_id -> [full LoFs, part LoFs, probability no LoF]
	HashMap<String, Double[]> gene_statistic;
	
	
	//sample_id -> fields of sample summary
	HashMap<String, SampleStat> sample_statistic;
	
	//transcript_id -> [LoFs, observed LoFs homo, observed LoFs hetero]
	HashMap<String, Integer[]> transcript_statistic;
	
	protected static final NodeLogger logger = NodeLogger.getLogger(LOFSummaryNodeModel.class);
	
	/**
	 * 
	 * @param vcf_file file containing annotation
	 * @param cds_file file containing gene id for each transcript id, can be found
	 *                 here: http://www.ensembl.org/info/data/ftp/index.html
	 */
	public Summarizer(String vcf_file, String cds_file) {
		this.vcf_file = vcf_file;
		this.cds_file = cds_file;
		try {
			logger.info("Reading CDS file...");
			this.readCDSFile();
			logger.info("CDS file read.");
		} catch (IOException e) {
			logger.error("Reading CDS file failed. Variant effect (full or partial) cannot be calculated.");
		}
		
		additional_titles = new ArrayList<String>();
		lof_statistic = new ArrayList<LoFVariant>();
	}
	
	public String getLoFStatistic() {
		if(lof_statistic.size()==0) {
			this.extract_LOFs();
		}
		String outfile = vcf_file.replace("vcf", "variant_summary.tsv");
		try {
			this.writeLOFStatistics(outfile);
			this.generateGeneStatistic();
			this.generateSampleStatistic();
		} catch (IOException e) {
			logger.error(outfile+ " could not be written!");
		}
		return outfile;
	}
	
	public String getGeneStatistic() {
		if(lof_statistic.size()==0) {
			this.getLoFStatistic();
		}
		String outfile = vcf_file.replace("vcf","gene_summary.tsv");
		try {
			this.writeGeneStatistic(outfile);
		} catch (IOException e) {
			logger.error(outfile+ " could not be written!");
		}
		return outfile;
	}
	
	public String getSampleStatistic() {
		if(lof_statistic.size()==0) {
			this.getLoFStatistic();
		}
		String outfile = vcf_file.replace("vcf", "sample_summary.tsv");
		try {
			this.writeSampleStatistic(outfile);
		} catch (IOException e) {
			logger.error(outfile+" could not be written!");
		}
		return outfile;
	}
	
//	public String getTranscriptStatistic() {
//		if(lof_statistic.size()==0) {
//			this.extract_LOFs();
//		}
//		String outfile = vcf_file.replace("vcf", "transcript_summary.tsv");
//		try {
//			this.generateTranscriptStatistic();
//			this.writeTranscriptStatistic(outfile);
//		} catch (IOException e) {
//			logger.error(outfile+" could not be written!");
//		}
//		return outfile;
//	}
	
//	public String getAnnotationStatistic() {
//		if(lof_statistic.size()==0) {
//			this.extract_LOFs();
//		}
//		String outfile = vcf_file.replace("vcf", "annotation_summary.tsv");
//		try {
//			this.writeAnnotationStats(outfile);
//		} catch (IOException e) {
//			logger.error(outfile+" could not be written!");
//		}
//		return outfile;
//	}
	
//	public String getKnockOutGenes() {
//		if(lof_statistic.size()==0) {
//			this.extract_LOFs();
//			this.generateGeneStatistic();
//			this.generateSampleStatistic();
//		}
//		String outfile = vcf_file.replace("vcf", "knockout_genes.tsv");
//		try {
//			this.writeKnockOutGenes(outfile);
//		} catch (IOException e) {
//			logger.error(outfile+" could not be written!");
//		}
//		return outfile;
//	}
	
//	public String getSampleKOGeneMatrix() {
//		if(lof_statistic.size()==0) {
//			this.extract_LOFs();
//			this.generateGeneStatistic();
//			this.generateSampleStatistic();
//		}
//		String outfile = vcf_file.replace("vcf", "sample_ko_genes.tsv");
//		try {
//			this.writeSampleKOGeneMatrix(outfile);
//		} catch (IOException e) {
//			logger.error(outfile+" could not be written!");
//		}
//		return outfile;
//	}
	
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
	
	abstract void getLoFVariant(String chr, String pos, String ref, String id, String alt, int GT_index,String[] infos, ArrayList<String> genotypes);
	
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
			String [] fields;
			String chr, pos, ref;
			String [] ids, alt_alleles, infos;
			String [] afs = null;
			ArrayList<String> genotypes;
			FileInputStream inputStream = null;
			Scanner sc = null;
			try {
			    inputStream = new FileInputStream(vcf_file);
			    sc = new Scanner(inputStream, "UTF-8");
			    while (sc.hasNextLine()) {
			        String line = sc.nextLine();
			        if(line.startsWith("#"))  {
			        	if(line.startsWith("##INFO=<ID=CSQ")) {
			        		vep_header = new HashMap<String,Integer>();
			        		String [] tmp = line.split("\\|");
			        		nrOfVEPfields = tmp.length;
			        		String part; 
			        		for(int i = 0; i<tmp.length; i++) {
			        			part = tmp[i];
			        			if(part.contains("ALLELE_NUM")) {
			        				vep_header.put("allele", i);
			        			} else if (part.contains("Consequence")) {
			        				vep_header.put("consequence", i);
			        			} else if (part.equals("SYMBOL")) {
			        				vep_header.put("gene_symbol",i);
			        			} else if (part.contains("Gene")) {
			        				vep_header.put("gene_id", i);
			        			} else if (part.contains("Feature")) {
			        				vep_header.put("transcript_id",i);
			        			} else if (part.contains("LoF_filter")) {
			        				additional_titles.add("lof_filter");
			        				vep_header.put("lof_filter",i);
			        			} else if (part.contains("LoF_flags")) {
			        				additional_titles.add("lof_flags");
			        				vep_header.put("lof_flags",i);
			        			} else if (part.contains("LoF_info")) {
			        				additional_titles.add("lof_info");
			        				vep_header.put("lof_info",i);
			        			} else if (part.contains("LoF")) {
			        				additional_titles.add("confidence");
			        				vep_header.put("confidence",i);
			        			} else if (part.equals("ExAC_AF")) {
			        				additional_titles.add("ExAC_AF");
			        				vep_header.put("ExAC_AF", i);
			        			} else if (part.equals("CADD_PHRED")) {
			        				additional_titles.add("CADD_PHRED");
			        				vep_header.put("CADD_PHRED",i);
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
			        	
			        	for(String i: infos) {
			        		if(i.startsWith("AF")) {
			        			afs = i.split("=")[1].split(",");
			        		}
			        	}
			        	
			        	for(int i = 0; i<alt_alleles.length;i++) {
			        		String id = ".";
			        		if(ids.length>i) {
			        			id = ids[i];
			        		}
			        		if(afs == null) {
			        			getLoFVariant(chr, pos, ref, id, alt_alleles[i], i+1,infos, genotypes);
			        		} else {
			        			getLoFVariant(chr, pos, ref, id, alt_alleles[i], i+1,infos, genotypes);
			        		}
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
			
		} catch (IOException e) {
			logger.error("VCF file couldn't be read!");
		}
	}
	
	private void generateGeneStatistic() {
		//gene_id -> sample_id -> [LoF transcripts homo, LoF transcripts hetereo]
		gene_sample_statistic = new HashMap<>();
		
		//gene_id -> sample_id -> is_unaffected
		gene_sample_unaffected = new HashMap<>();
		for(String gene: gene_id2gene_symbol.keySet()) {
			for(String sample_id : sample_ids ) {
				HashMap<String,Boolean> my = gene_sample_unaffected.get(gene);
				if(my == null) my = new HashMap<>();
				my.put(sample_id,true);
				gene_sample_unaffected.put(gene, my);
			}
		}
		
		//gene_id -> sample_id -> transcripts
		HashMap<String, HashMap<String, HashSet<String>>> transcripts_hom = new HashMap<>();
		HashMap<String, HashMap<String, HashSet<String>>> transcripts_het = new HashMap<>();
		
		//temporary transcript sets for phased genotypes
		HashMap<String, HashMap<String, HashSet<String>>> transcripts_pat = new HashMap<>();
		HashMap<String, HashMap<String, HashSet<String>>> transcripts_mat = new HashMap<>();
		
		//gene_id -> [full LoFs, part LoFs]
		gene_statistic = new HashMap<>();
		
		for(LoFVariant lof: lof_statistic)  {
			HashMap<String, LoFGene> gene2lof = lof.getGene_id2lofgene_fields();
			
			for(String gene_id: gene2lof.keySet()) {
				//initialize gene statistic
				Double [] stats = gene_statistic.get(gene_id);
				if(stats==null) {
					stats = new Double[3];
					for(int i=0; i<2;i++) {
						stats[i]=0.0;
					}
					stats[2] = 1.0;
				}
				
				//calculate probability that gene contains no LoF
				stats[2] = stats[2] * Math.pow((1.0 - lof.getAF()),2);
				
				ArrayList<String> transList = gene2lof.get(gene_id).getTranscripts();
				//copy to remove duplicates
				HashSet<String> transcripts = new HashSet<>();
				for(String t:transList) {
					transcripts.add(t);
				}
				//is variant full LOF
				if(transcripts.size() == gene_id2transcript_ids.get(gene_id).size()) {
					stats[0]++;
				} else {
					stats[1]++;
				}
				
				for(int i=0; i<lof.getGenotypes().size(); i++) {
					String gt = lof.getGenotypes().get(i);
					
					if(gt.contains(".") || (gt.charAt(0)>'0' || gt.charAt(2)>'0')) {
						HashMap<String,Boolean> my = gene_sample_unaffected.get(gene_id);
						if(my == null) my = new HashMap<>();
						my.put(sample_ids.get(i), false);
						gene_sample_unaffected.put(gene_id, my);
					}
					
					//homozygous affected transcripts
					if(gt.charAt(0)==gt.charAt(2) && gt.charAt(0)>'0') {
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
					} else if(gt.contains("|")) {
						//phased heterozygous phenotype
						if(gt.charAt(0)>'0') {
							//add to transcripts_mat
							HashMap<String, HashSet<String>> sample2transcripts = transcripts_mat.get(gene_id);
							if(sample2transcripts==null) {
								sample2transcripts = new HashMap<>();
							}
							HashSet<String> tmp_set = sample2transcripts.get(sample_ids.get(i));
							if(tmp_set == null) {
								tmp_set = new HashSet<>();
							}
							
							tmp_set.addAll(transcripts);
							sample2transcripts.put(sample_ids.get(i), tmp_set);
							transcripts_mat.put(gene_id, sample2transcripts);
						} else if (gt.charAt(2)>'0') {
							//add to transcripts_pat
							HashMap<String, HashSet<String>> sample2transcripts = transcripts_pat.get(gene_id);
							if(sample2transcripts==null) {
								sample2transcripts = new HashMap<>();
							}
							HashSet<String> tmp_set = sample2transcripts.get(sample_ids.get(i));
							if(tmp_set == null) {
								tmp_set = new HashSet<>();
							}
							
							tmp_set.addAll(transcripts);
							sample2transcripts.put(sample_ids.get(i), tmp_set);
							transcripts_pat.put(gene_id, sample2transcripts);
						}
					} else if (gt.charAt(0)>'0' || gt.charAt(2)>'0') {
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
		
		//merge transcripts_mat & transcripts_pat
		
		int hom_lof_nr, het_lof_nr;
		for(String gene: gene_id2gene_symbol.keySet()) {
			for(String sample: sample_ids) {
				HashSet<String> hom_lofs = new HashSet<>();
				HashSet<String> het_lofs = new HashSet<>();
				
				if(transcripts_hom.get(gene)!=null) {
					hom_lofs = transcripts_hom.get(gene).get(sample);
				}
				if(transcripts_het.get(gene)!=null) {
					het_lofs = transcripts_het.get(gene).get(sample);
				} else if(transcripts_mat.get(gene)!=null || transcripts_pat.get(gene)!=null) {
					
					if(transcripts_mat.get(gene)!=null && transcripts_pat.get(gene)!=null) {
						//overlap is compound heterozygous
						HashSet<String> tmp_set_mat, tmp_set_pat;
						tmp_set_mat = transcripts_mat.get(gene).get(sample);
						tmp_set_pat = transcripts_pat.get(gene).get(sample);
						
						if(tmp_set_mat==null) tmp_set_mat = new HashSet<>();
						if(tmp_set_pat==null) tmp_set_pat = new HashSet<>();
						
						HashSet<String> intersection = new HashSet<>(tmp_set_mat);

						intersection.retainAll(tmp_set_pat);
						//TESTING
//						if(intersection.size()>0) {
//							System.out.println("Compound heterozygotes detected for gene: "+gene+" in sample: "+sample);
//							System.out.println("Affected transcripts");
//							for(String i:intersection) {
//								System.out.println(i);
//							}
//						}
						if(hom_lofs == null) hom_lofs = new HashSet<>();
						hom_lofs.addAll(intersection);
						
						HashSet<String> differenceA = new HashSet<>(tmp_set_mat);
						differenceA.removeAll(tmp_set_pat);
						het_lofs.addAll(differenceA);
						
						HashSet<String> differenceB = new HashSet<>(tmp_set_pat);
						differenceB.removeAll(tmp_set_mat);
						het_lofs.addAll(differenceB);
					} else if (transcripts_mat.get(gene)!=null) {
						HashSet<String> tmp = transcripts_mat.get(gene).get(sample);
						if(tmp == null) tmp = new HashSet<>();
						het_lofs.addAll(tmp);
					} else if (transcripts_pat.get(gene)!=null){
						HashSet<String> tmp = transcripts_pat.get(gene).get(sample);
						if(tmp == null) tmp = new HashSet<>();
						het_lofs.addAll(tmp);
					}
				}
				
				if(het_lofs==null) het_lofs = new HashSet<>();
				if(hom_lofs==null) hom_lofs = new HashSet<>();
				
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
		//sample id -> number of LoF variants
		HashMap<String, Integer> fullLOFs = new HashMap<>();
		HashMap<String, Integer> partLOFs = new HashMap<>();
		boolean isFull = false;
		
		//initialize fullLOFs and partLOFs
		for(String sample:sample_ids) {
			fullLOFs.put(sample, 0);
			partLOFs.put(sample, 0);
		}
		
		//fill fullLOFs and partLOFs
		
		//iterate over all LoF variants
		for(LoFVariant lof:lof_statistic) {
			HashMap<String, LoFGene> lof_genes = lof.getGene_id2lofgene_fields();
			isFull = false;
			
			//check whether variant effect is full or partial
			for(Entry<String, LoFGene> gene:lof_genes.entrySet()) {
				int all_trans = gene_id2transcript_ids.get(gene.getKey()).size();
				int lof_trans = gene.getValue().getNrOfTranscripts();
				if(all_trans == lof_trans) {
					isFull = true;
				}
			}
			
			//iterate over all samples/genotypes
			for(int i = 0; i< lof.getGenotypes().size(); i++) {
				String gt = lof.getGenotypes().get(i);
				String sample = sample_ids.get(i);
				
				//check whether sample is affected
				if(gt.charAt(0)>'0' || gt.charAt(2)>'0') {
					if(isFull) {
						Integer t = fullLOFs.get(sample);
						t++;
						fullLOFs.put(sample, t);
					} else {
						Integer x = partLOFs.get(sample);
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
				int lof_het = gene_sample_statistic.get(gene).get(sample)[1];

				if(all_trans == lof_hom) {//all transcripts are affected by LoF variants
					completeLOFgene.add(gene); 
				} else if((lof_hom>0) ||(lof_het>0)){
					partLOFgene.add(gene);
				}
			}
			sample_statistic.put(sample, new SampleStat(fullLOFs.get(sample),partLOFs.get(sample), completeLOFgene, partLOFgene));
		}
	}
	
//	private void generateTranscriptStatistic() {
//		transcript_statistic = new HashMap<>();
//		for(LoFVariant lofvar: lof_statistic) {
//			HashSet<String> transcripts = new HashSet<>();
//			HashMap<String, LoFGene> lofgene = lofvar.getGene_id2lofgene_fields();
//			int obs_hom = 0;
//			int obs_het = 0;
//			for(String g: lofgene.keySet()) {
//				transcripts.addAll(lofgene.get(g).getTranscripts());
//			}
//			for(String gt: lofvar.getGenotypes()) {
//				if(gt.charAt(0)>'0' && gt.charAt(2)==gt.charAt(0)) {
//					obs_hom ++;
//				} else if(gt.charAt(0)>'0' || gt.charAt(2)>'0') {
//					obs_het ++;
//				}
//			}
//			for(String t: transcripts) {
//				Integer[] stats = transcript_statistic.get(t);
//				if(stats==null) {
//					stats = new Integer[3];
//					for (int i = 0; i< 3; i++) {
//						stats[i]= 0;
//					}
//				}
//				stats[0]++;
//				stats[1]+= obs_hom;
//				stats[2]+= obs_het;
//				transcript_statistic.put(t, stats);
//			}
//		}
//	}
	
	private void writeLOFStatistics(String outfile) throws IOException{
		
		BufferedWriter bw = new BufferedWriter(new FileWriter(outfile));
		String [] header = {"chr", "pos", "rsId", "ref_allele", "alt_allele", "AF", "gene_id", "gene_symbol", "effect", "consequence", "lof_trans", "all_trans", "obs_hom", "obs_het"};
		
		//write header
		bw.write(header[0]);
		for(int i = 1; i < header.length; i++) {
			bw.write("\t"+header[i]);
		}
		for(String s: additional_titles) {
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
					logger.error("According to the CDS file gene "+s+" has no protein coding transcripts. The reference genome version of the CDS file and of the files used for the annotation have to be the equal.");
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
					ArrayList<String> infos = lofgene.getInfo(i);
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
			
			bw.write(l.getChr()+"\t"+l.getPos()+"\t"+l.getRsId()+"\t"+l.getRef_allele()+"\t"+l.getAlt_allele()+"\t"+l.getAF()+"\t"+gene_ids+"\t"+genes+"\t"+effects+"\t"+consequences+"\t"+lof_trans+"\t"+all_trans+"\t"+l.getObserved_homo()+"\t"+l.getObserved_hetero()+"\t");
			
			//write additional fields
			for(String s: add_infos) {
				s = s.replaceFirst(",","");
				bw.write(s+"\t");
			}
			bw.newLine();
		}
		
		bw.close();
	}
	
	private void writeGeneStatistic(String outfile) throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(outfile));
//		bw.write("#analyzed samples="+sample_ids.size());
//		bw.newLine();
		bw.write("gene_id\tgene_symbol\tfull\tpartial\tP(LoF>=1)\tunaffected_samples\taffected_samples\tko_samples\tobserved_LoF_frequency");
		bw.newLine();
		for(String gene: gene_id2gene_symbol.keySet()) {
			bw.write(gene+"\t"+gene_id2gene_symbol.get(gene));
			
			for(int i = 0; i<2 ; i++) {
				bw.write("\t"+gene_statistic.get(gene)[i].intValue());
			}
			
			bw.write("\t"+(1-gene_statistic.get(gene)[2]));
			
			HashSet<String> affected_samples = new HashSet<>();
			HashSet<String> ko_samples = new HashSet<>();
			int healthy_samples = 0;
			
			for(String sample_id: sample_ids) {
				if(gene_sample_unaffected.get(gene).get(sample_id)) {
					healthy_samples++;
				}
				if(sample_statistic.get(sample_id).getPart_LOF_genes().contains(gene)) {
					affected_samples.add(sample_id);
				} else if (sample_statistic.get(sample_id).getComplete_LOF_genes().contains(gene)) {
					ko_samples.add(sample_id);
				}
			}
			bw.write("\t"+healthy_samples);
			int affected = affected_samples.size()+ko_samples.size();
			bw.write("\t"+affected);
//			String sample_ids = "";
//			for(String id: affected_samples) {
//				sample_ids += ","+id;
//			}
//			sample_ids = sample_ids.replaceFirst(",", "");
//			bw.write("\t"+sample_ids);
			
			bw.write("\t"+ko_samples.size());
//			sample_ids = "";
//			for(String id: ko_samples) {
//				sample_ids += ","+id;
//			}
//			sample_ids = sample_ids.replaceFirst(",", "");
//			bw.write("\t"+sample_ids);
			
			double frequency = ((double) affected)/((double) (healthy_samples+affected));
			bw.write("\t"+frequency);
			bw.newLine();
		}
		
		bw.close();
	}
	
	private void writeSampleStatistic(String outfile) throws IOException {
		
		BufferedWriter bw = new BufferedWriter(new FileWriter(outfile));
		bw.write("sample_id\tfull\tpartial\taffectedGenes\tknocked_out\tcompleteLOFgenes\taffectedLOFGenes");
		bw.newLine();
		for(String s: sample_statistic.keySet()) {
			SampleStat stat = sample_statistic.get(s);
			ArrayList<String> completes = stat.getComplete_LOF_genes();
			Collections.sort(completes);
			int affected = stat.getPart_LOF_genes().size()+completes.size();
			bw.write(s+"\t"+stat.getFullLOFs()+"\t"+stat.getPartLOFs()+"\t"+affected+"\t"+completes.size()+"\t");
			
			String affectedGeneList = "";
			
			if(completes.size() >= 1) {
				String gene = completes.get(0);
				affectedGeneList = gene;
				bw.write(gene +":"+gene_id2gene_symbol.get(gene));
			}
			
			for(int i = 1; i< completes.size(); i++) {
				String u = completes.get(i);
				affectedGeneList += ","+u;
				bw.write(","+u+":"+gene_id2gene_symbol.get(u));
			}
			
			for(String affGene: stat.getPart_LOF_genes()) {
				if(affectedGeneList.equals("")) {
					affectedGeneList += affGene;
				}
				affectedGeneList += ","+affGene;
			}
			
			bw.write("\t"+affectedGeneList);
			bw.newLine();
		}
		bw.close();
		//end of original method
	}
	
//	private void writeKnockOutGenes (String outfile) throws IOException {
//		HashSet<String> knockout_genes = new HashSet<>();
//		for(String s: sample_statistic.keySet()) {
//			knockout_genes.addAll(sample_statistic.get(s).getComplete_LOF_genes());
//		}
//		
//		BufferedWriter bw = new BufferedWriter(new FileWriter(outfile));
//		for(String s: knockout_genes) {
//			bw.write(s+"\t"+this.gene_id2gene_symbol.get(s));
//			bw.newLine();
//		}
//		bw.close();
//	}
	
//	private void writeSampleKOGeneMatrix(String outfile) throws IOException {
//		HashSet<String> knockout_genes = new HashSet<>();
//		for(String s: sample_statistic.keySet()) {
//			knockout_genes.addAll(sample_statistic.get(s).getComplete_LOF_genes());
//		}
//		
//		Object [] all_ko_genes = knockout_genes.toArray();
//		
//		String header = "";
//		for(Object gene: all_ko_genes) {
//			header += "\t" + this.gene_id2gene_symbol.get(gene+"");
//		}
//		header = header.replaceFirst("\t", "");
//		
//		String line = "";
//		
//		BufferedWriter bw = new BufferedWriter(new FileWriter(outfile));
//		bw.write(header);
//		bw.newLine();
//		for(String s: sample_ids) {
//			line = s;
//			ArrayList<String> ko_genes = sample_statistic.get(s).getComplete_LOF_genes();
//			for(Object gene: all_ko_genes) {
//				if(ko_genes.contains(gene+"")) {
//					line += "\t1";
//				} else {
//					line += "\t0";
//				}
//			}
//			bw.write(line);
//			bw.newLine();
//		}
//		
//		
//		bw.close();
//	}
	
//	private void writeAnnotationStats (String outfile) throws IOException {
//		int hc_lof=0;
//		int all_lof=0;
//		TreeMap<String, Integer> term2hc_lof = new TreeMap<>();
//		TreeMap<String, Integer> term2all_lof = new TreeMap<>();
//		
//		int conf_index = additional_titles.indexOf("confidence");
//		
//		for(LoFVariant lof: lof_statistic) {
//			for(LoFGene gene: lof.getGene_id2lofgene_fields().values()) {
//				ArrayList<String> consequences = gene.getConsequences();
//				for(int i = 0; i < consequences.size(); i++) {
//					String term = consequences.get(i);
//					Integer count = term2all_lof.get(term);
//					if(count == null) {
//						count = 0;
//					}
//					count ++;
//					term2all_lof.put(term, count);
//					all_lof++;
//					if(conf_index != -1 && gene.getInfo(conf_index).get(i).equals("HC")) {
////						System.out.println(lof.getChr()+" "+lof.getPos()+" "+lof.getRef_allele()+" "+term+" HC gefunden");
//						hc_lof++;
//						Integer hc_count = term2hc_lof.get(term);
//						if(hc_count == null) {
//							hc_count = 0;
//						}
//						hc_count++;
//						term2hc_lof.put(term, hc_count);
//					}
//				}
//			}
//		}
//		
//		BufferedWriter bw = new BufferedWriter(new FileWriter(outfile));
//		bw.write("Number of LoF annotations:\t"+all_lof);
//		bw.newLine();
//		bw.write("Number of hc LoF annotations:\t"+hc_lof);
//		bw.newLine();
//		
//		bw.write("#Number of LoF annotations per term:");
//		bw.newLine();
//		for(Map.Entry<String, Integer> entry: term2all_lof.entrySet()) {
//			bw.write(entry.getKey()+"\t"+entry.getValue());
//			bw.newLine();
//		}
//		
//		bw.write("#Number of hc LoF annotations per term:");
//		bw.newLine();
//		for(Map.Entry<String, Integer> entry: term2hc_lof.entrySet()) {
//			bw.write(entry.getKey()+"\t"+entry.getValue());
//			bw.newLine();
//		}
//		bw.close();
//		
//	}
	
//	private void writeTranscriptStatistic(String outfile) throws IOException {
//		BufferedWriter bw = new BufferedWriter(new FileWriter(outfile));
//		
//		bw.write("transcript_id\tgene_id\tgene_symbol\tlofs\tobs_hom\tobs_het");
//		bw.newLine();
//		
//		for(String t: transcript_statistic.keySet()) {
//			String gene_id = transcript_id2gene_id.get(t);
//			String gene_symbol = gene_id2gene_symbol.get(gene_id);
//			Integer [] stats = transcript_statistic.get(t);
//			bw.write(t+"\t"+gene_id+"\t"+gene_symbol+"\t"+stats[0]+"\t"+stats[1]+"\t"+stats[2]);
//			bw.newLine();
//		}
//		
//		bw.close();
//	}
	
}