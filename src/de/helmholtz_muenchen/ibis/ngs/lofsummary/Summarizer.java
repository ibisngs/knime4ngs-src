package de.helmholtz_muenchen.ibis.ngs.lofsummary;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
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
	String ped_file;
	
	//summaries
	
	//variant summary
	ArrayList<LoFVariant> lof_statistic;
	ArrayList<String> additional_titles;
	HashMap<String,Integer> vep_header; //contains indices for LoF, LoF_filter etc.
	int nrOfVEPfields;
	
	//gene_id+sample_id as key
	HashMap<String, GeneSampleInfo> gene_sample_statistic;
	
	HashMap<String, GeneInfo> gene_statistic;
	HashMap<String, HashSet<String>> gene2transcripts;
	
	//save order of sample_ids
	ArrayList<String> my_sample_ids;
	HashMap<String, SampleInfo> sample_statistic;
	
	protected static final NodeLogger logger = NodeLogger.getLogger(LOFSummaryNodeModel.class);
	
	/**
	 * 
	 * @param vcf_file file containing annotation
	 * @param cds_file file containing gene id for each transcript id, can be found
	 *                 here: http://www.ensembl.org/info/data/ftp/index.html
	 */
	public Summarizer(String vcf_file, String cds_file, String ped_file) {
		this.vcf_file = vcf_file;
		this.cds_file = cds_file;
		this.ped_file = ped_file;
		sample_statistic = new HashMap<>();
		gene_statistic = new HashMap<>();
		gene_sample_statistic = new HashMap<>();
		gene2transcripts = new HashMap<>();
		try {
			logger.info("Reading CDS file...");
			this.readCDSFile();
			this.readPEDFile();
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
			this.generateGeneSampleStatistic();
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
	
	public String getTrioStatistic() {
		if(lof_statistic.size()==0) {
			this.getLoFStatistic();
		}
		String outfile = vcf_file.replace("vcf", "trio_summary.tsv");
		try {
			this.writeTrioSummary(outfile);
		} catch (IOException e) {
			logger.error(outfile+" could not be written!");
		}
		return outfile;
	}
	
	/**
	 * reads CDS file and fills transcript_id2gene_id and gene_id2transcript_ids
	 * @throws IOException
	 */
	private void readCDSFile() throws IOException {
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
		        	
		        	if(gene2transcripts.containsKey(gene_id)) {
		        		gene2transcripts.get(gene_id).add(transcript_id);
		        	} else {
		        		HashSet<String> tmp = new HashSet<>();
		        		tmp.add(transcript_id);
		        		gene2transcripts.put(gene_id,tmp);
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
	
	private void readPEDFile() {
		
		
		try {
			BufferedReader br = new BufferedReader(new FileReader(this.ped_file));
			String line;
			String [] fields;
			while((line=br.readLine())!=null) {
				if(line.startsWith("Family ID")) continue;
				
				fields = line.split("\t");
				if(!sample_statistic.containsKey(fields[1])) {
					SampleInfo sample_info = new SampleInfo(fields[0],fields[2],fields[3],fields[4],fields[5]);
					sample_statistic.put(fields[1],sample_info);
				}
			}
		} catch (FileNotFoundException e) {
			logger.error(ped_file+" not found!");
			e.printStackTrace();
		} catch (IOException e) {
			logger.error(ped_file+" could not be read!");
			e.printStackTrace();
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
			        		my_sample_ids = new ArrayList<>();
			        		for(int i = 9; i< header.length; i++) {
			        			my_sample_ids.add(header[i]);
			        			if(!sample_statistic.containsKey(header[i])) {
			        				sample_statistic.put(header[i], new SampleInfo());
			        			}
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
	
	private void generateGeneSampleStatistic() {
		
		for(LoFVariant lof: lof_statistic)  {
			HashMap<String, LoFGene> gene2lof = lof.getGene_id2lofgene_fields();
			
			for(String gene_id: gene2lof.keySet()) {
				
				
				GeneInfo gene_info = gene_statistic.get(gene_id);
				
				gene_info.addProb(lof.getPos(), lof.getAF());
				
				ArrayList<String> transList = gene2lof.get(gene_id).getTranscripts();
				//copy to remove duplicates
				HashSet<String> transcripts = new HashSet<>();
				for(String t:transList) {
					transcripts.add(t);
				}
				//is variant full LOF
				if(transcripts.size() == gene2transcripts.get(gene_id).size()) {
					gene_info.incrementFullLoFs();
				} else {
					gene_info.incrementPartLoFs();
				}
				
				for(int i=0; i<lof.getGenotypes().size(); i++) {
					String gt = lof.getGenotypes().get(i);
					
					GeneSampleInfo gsi =  gene_sample_statistic.get(gene_id+"_"+my_sample_ids.get(i));
					if(gsi == null) {
						gsi = new GeneSampleInfo();
						gene_sample_statistic.put(gene_id+"_"+my_sample_ids.get(i), gsi);
					}
					
					if(gt.contains(".") || (gt.charAt(0)>'0' || gt.charAt(2)>'0')) {
						gsi.setAffected();
					}
					
					//homozygous affected transcripts
					if(gt.charAt(0)==gt.charAt(2) && gt.charAt(0)>'0') {
						gsi.addHomTranscripts(transcripts);
					} else if(gt.contains("|")) {
						//phased heterozygous phenotype
						if(gt.charAt(0)>'0') {
							gsi.addMatTranscripts(transcripts);
						} else if (gt.charAt(2)>'0') {
							gsi.addPatTranscripts(transcripts);
						}
					} else if (gt.charAt(0)>'0' || gt.charAt(2)>'0') {
						gsi.addHetTranscripts(transcripts);
					}	
				}
			}
		}
		
		for(String key : gene_sample_statistic.keySet()) {
				gene_sample_statistic.get(key).calculate();
		}
	}
	
	private void generateSampleStatistic() {

		boolean isFull = false;
		
		//fill fullLOFs and partLOFs
		
		//iterate over all LoF variants
		for(LoFVariant lof:lof_statistic) {
			HashMap<String, LoFGene> lof_genes = lof.getGene_id2lofgene_fields();
			isFull = false;
			//check whether variant effect is full or partial
			for(Entry<String, LoFGene> gene:lof_genes.entrySet()) {
				int all_trans = gene2transcripts.get(gene.getKey()).size();
				int lof_trans = gene.getValue().getNrOfTranscripts();
				if(all_trans == lof_trans) {
					isFull = true;
				}
			}
			
			//iterate over all samples/genotypes
			for(int i = 0; i< lof.getGenotypes().size(); i++) {
				String gt = lof.getGenotypes().get(i);
				String sample = my_sample_ids.get(i);
				
				//check whether sample is affected
				if(gt.charAt(0)>'0' || gt.charAt(2)>'0') {
					if(isFull) {
						sample_statistic.get(sample).incrementFullLOFs();
					} else {
						sample_statistic.get(sample).incrementPartLOFs();
					}
				}
			}
		}
		
		//get complete LOF genes
		for(Entry<String,GeneSampleInfo> e: gene_sample_statistic.entrySet()) {
			int all_trans = gene2transcripts.get(e.getKey().split("_")[0]).size();
			int lof_hom = e.getValue().getHomTrans();
			int lof_het = e.getValue().getHetTrans();
			
			SampleInfo si = sample_statistic.get(e.getKey().split("_")[1]);
			if(all_trans == lof_hom) {//all transcripts are affected by LoF variants
				si.addCompleteLoFGene(e.getKey().split("_")[0]);
			} else if((lof_hom>0) ||(lof_het>0)){
				si.addPartLoFGene(e.getKey().split("_")[0]);
			}
		}
	}
	
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
				genes += ","+ gene_statistic.get(s).getSymbol();
				String nrOfTrans = "-";
				if(gene2transcripts.containsKey(s)) {
					nrOfTrans = gene2transcripts.get(s).size()+"";	
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
		bw.write("gene_id\tgene_symbol\tfull\tpartial\tP(LoF>=1)\tunaffected_samples\taffected_samples\tko_samples\tobserved_LoF_frequency");
		bw.newLine();
		for(String gene: gene_statistic.keySet()) {

			GeneInfo gi = gene_statistic.get(gene);
			bw.write(gene+"\t"+gi.getSymbol());
			bw.write("\t"+gi.getFullLoFs());
			bw.write("\t"+gi.getPartLoFs());
			bw.write("\t"+(1-gi.getProbUnaffected()));
			
			HashSet<String> affected_samples = new HashSet<>();
			HashSet<String> ko_samples = new HashSet<>();
			int healthy_samples = 0;
			
			for(String sample_id: my_sample_ids) {
				if(gene_sample_statistic.get(gene+"_"+sample_id).isUnaffected()) {
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
			bw.write("\t"+ko_samples.size());
			
			double frequency = ((double) affected)/((double) (healthy_samples+affected));
			bw.write("\t"+frequency);
			bw.newLine();
		}
		
		bw.close();
	}
	
	private void writeSampleStatistic(String outfile) throws IOException {
		
		HashSet<String> ko_genes = new HashSet<>();
		
		BufferedWriter bw = new BufferedWriter(new FileWriter(outfile));
		bw.write("sample_id\tfull\tpartial\taffectedGenes\tknocked_out\tcompleteLOFgenes\taffectedLOFGenes");
		bw.newLine();
		for(String s: sample_statistic.keySet()) {
			SampleInfo stat = sample_statistic.get(s);
			ArrayList<String> completes = stat.getComplete_LOF_genes();
			Collections.sort(completes);
			int affected = stat.getPart_LOF_genes().size()+completes.size();
			bw.write(s+"\t"+stat.getFullLOFs()+"\t"+stat.getPartLOFs()+"\t"+affected+"\t"+completes.size()+"\t");
			
			String affectedGeneList = "";
			
			if(completes.size() >= 1) {
				String gene = completes.get(0);
				ko_genes.add(gene);
				affectedGeneList = gene;
				bw.write(gene +":"+gene_statistic.get(gene).getSymbol());
			}
			
			for(int i = 1; i< completes.size(); i++) {
				String u = completes.get(i);
				ko_genes.add(u);
				affectedGeneList += ","+u;
				bw.write(","+u+":"+gene_statistic.get(u).getSymbol());
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
		writeGeneList(ko_genes,outfile.replace("sample_summary.tsv", "ko_genes.tsv"));
	}
	
	private void writeTrioSummary(String outfile) throws IOException {
		
		String header = "sample_id\tde_novo_LOF\tKOs_de_novo\tde_novo_KO";
		ArrayList<String> lof_genes = new ArrayList<>();
		ArrayList<String> ko_genes = new ArrayList<>();
		
		HashSet<String> de_novo_LOF_genes = new HashSet<>();
		HashSet<String> de_novo_KO_genes = new HashSet<>();
		
		BufferedWriter bw = new BufferedWriter(new FileWriter(outfile));
		bw.write(header);
		bw.newLine();
		for(String s:sample_statistic.keySet()) {
			SampleInfo stat = sample_statistic.get(s);
			lof_genes = stat.getPart_LOF_genes();
			lof_genes.addAll(stat.getComplete_LOF_genes());
			ko_genes = stat.getComplete_LOF_genes();
			boolean has_parents = false;
			
			if(my_sample_ids.contains(stat.getMatId())) {
				has_parents = true;
				SampleInfo mat = sample_statistic.get(stat.getMatId());
				lof_genes.removeAll(mat.complete_LOF_genes);
				lof_genes.removeAll(mat.part_LOF_genes);
				
				ko_genes.removeAll(mat.complete_LOF_genes);
			}
			
			if(my_sample_ids.contains(stat.getPatId())) {
				has_parents = true;
				SampleInfo pat = sample_statistic.get(stat.getPatId());
				lof_genes.removeAll(pat.complete_LOF_genes);
				lof_genes.removeAll(pat.part_LOF_genes);
				
				ko_genes.removeAll(pat.complete_LOF_genes);
			}
			
			if(has_parents) {
			
				String de_novo_LOF = "";
				for(String d: lof_genes) {
					de_novo_LOF_genes.add(d);
					de_novo_LOF += ","+d;
				}
				de_novo_LOF = de_novo_LOF.replaceFirst(",", "");
			
				String de_novo_KO = "";
				for(String d: ko_genes) {
					de_novo_KO_genes.add(d);
					de_novo_KO += ","+d+":"+gene_statistic.get(d).getSymbol();
				}
				de_novo_KO = de_novo_KO.replaceFirst(",","");
			
				bw.write(s);
				bw.write("\t"+lof_genes.size());
				bw.write("\t"+ko_genes.size());
				bw.write("\t"+de_novo_KO);
				bw.newLine();
			}
		}
		bw.close();
		writeGeneList(de_novo_LOF_genes,outfile.replace("trio_summary.tsv","de_novo_LOF_genes.tsv"));
		writeGeneList(de_novo_KO_genes,outfile.replace("trio_summary.tsv","de_novo_KO_genes.tsv"));
	}
	
	private void writeGeneList(HashSet<String> genes, String outfile) throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(outfile));
		for(String gene: genes) {
			
			bw.write(gene+"\t"+gene_statistic.get(gene).getSymbol());
			
			String samples = "";
			for(String sample_id : my_sample_ids) {
				if(sample_statistic.get(sample_id).getPart_LOF_genes().contains(gene) || sample_statistic.get(sample_id).getComplete_LOF_genes().contains(gene)) {
					samples += ","+sample_id;
				}
			}
			
			bw.write("\t"+samples.replaceFirst(",", ""));
			bw.newLine();
		}
		bw.close();
	}
}