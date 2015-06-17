package de.helmholtz_muenchen.ibis.ngs.lofsummary;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;
import java.util.Scanner;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.knime.core.node.NodeLogger;

public abstract class Summarizer {
	
	String vcf_file;
	String cds_file;
	String ped_file;
	String geneback_file;
	ArrayList<String> geneback_header= new ArrayList<>();
	
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
	
	//gene_id -> [dev_aff, dev_in, called_samples]
	HashMap<String, Double[]> gene2background;
	
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
	public Summarizer(String vcf_file, String cds_file, String ped_file, String geneback_file) {
		
		this.vcf_file = vcf_file;
		this.cds_file = cds_file;
		this.ped_file = ped_file;
		this.geneback_file = geneback_file;
		
		sample_statistic = new HashMap<>();
		gene_statistic = new HashMap<>();
		gene_sample_statistic = new HashMap<>();
		gene2transcripts = new HashMap<>();
		gene2background = new HashMap<>();
		additional_titles = new ArrayList<String>();
		lof_statistic = new ArrayList<LoFVariant>();
		
		if(cds_file != null) {
			try {
				logger.info("Reading CDS file...");
				this.readCDSFile();
				logger.info("CDS file read.");
			} catch (IOException e) {
				logger.warn("Reading CDS file failed. Variant effect (full or partial) cannot be calculated.");
			}
		}
		
		if(ped_file != null) {
			try {
				logger.info("Reading PED file...");
				this.readPEDFile();
				logger.info("PED file read.");
			} catch (IOException e) {
				logger.warn("Reading PED file failed. De novo LoF genes and KO genes cannot be found.");
			}
		}
		
		if(geneback_file != null) {
			try {
				logger.info("Reading genetic background file...");
				this.readGeneBackFile();
				logger.info("Genetic background file read.");
			} catch (IOException e) {
				logger.warn("Reading genetic background file failed.");
			}
		}
	}
	
	public String [] getSummaries() {

		String [] result = new String[5];
		try {
			logger.info("Reading VCF file...");
			this.extract_LOFs();
			
			logger.info("Writing variant summary...");
			String outfile1 = vcf_file.replace("vcf", "variant_summary.tsv");
			result[0] = outfile1;
			this.writeLOFStatistics(outfile1);
			
			this.generateSampleStatistic();
			
			logger.info("Generating gene and sample summaries...");
			this.generateGeneSampleStatistic();

			logger.info("Writing gene summary...");
			String outfile2 = vcf_file.replace("vcf","gene_summary.tsv");
			result[1] = outfile2;
			this.writeGeneStatistic(outfile2);
			
			logger.info("Writing sample summary...");
			String outfile3 = vcf_file.replace("vcf", "sample_summary.tsv");
			result[2] = outfile3;
			this.writeSampleStatistic(outfile3);
			
			String outfile4 = "";
			if(ped_file!=null) {
				logger.info("Writing trio summary...");
				outfile4 = vcf_file.replace("vcf", "trio_summary.tsv");
				this.writeTrioSummary(outfile4);
			}
			result[3]=outfile4;
			
			logger.info("Writing genetic background...");
			String outfile5 = vcf_file.replace("vcf", "gene_back.tsv");
			result[4] = outfile5;
			this.writeGeneticBackground(outfile5);
			
		} catch (IOException e) {
			logger.error(e.getMessage());
		}
		return result;
	}
	
	/**
	 * reads CDS file and fills gene_id2transcript_ids
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
	
	private void readPEDFile() throws IOException{
		
		BufferedReader br = new BufferedReader(new FileReader(this.ped_file));
		String line;
		String[] fields;
		while ((line = br.readLine()) != null) {
			if (line.startsWith("Family ID"))
				continue;

			fields = line.split("\t");
			if (!sample_statistic.containsKey(fields[1])) {
				SampleInfo sample_info = new SampleInfo(fields[0], fields[2],fields[3], fields[4], fields[5]);
				sample_statistic.put(fields[1], sample_info);
			}
		}
		br.close();
	}
	
	private void readGeneBackFile() throws IOException {
		
		BufferedReader br = new BufferedReader(new FileReader(this.geneback_file));
		String line;
		String[] fields;
		while((line = br.readLine())!=null) {
			if(line.startsWith("#")) {
				geneback_header.add(line);
				continue;
			}
			fields = line.split("\t");
			Double[] counts = new Double[3];
			for(int i = 0; i < 3; i++) {
				counts[i] = Double.parseDouble(fields[i+1]);
			}
			gene2background.put(fields[0], counts);
			
		}
		br.close();
	}
	
	abstract void getLoFVariant(String chr, String pos, String ref, String id, String alt, double af, int GT_index,String[] infos, ArrayList<String> genotypes);
	
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
			String[] ac_adj = null;
			String[] ac = null;
			String an_adj = null;
			String an = null;
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
			        	
			        	boolean adj_not_found = true;
			        	for(String i: infos) {
			        		if(i.startsWith("AC_Adj=")) {
			        			adj_not_found = false;
			        			ac_adj = i.split("=")[1].split(",");
			        		} else if (i.startsWith("AN_Adj=")) {
			        			an_adj = i.split("=")[1];
			        		} else if (i.startsWith("AN=")) {
			        			an = i.split("=")[1];
			        		} else if (i.startsWith("AC=")) {
			        			ac = i.split("=")[1].split(",");
			        		}
			        	}
			        	
			        	for(int i = 0; i<alt_alleles.length;i++) {
			        		String id = ".";
			        		if(ids.length>i) {
			        			id = ids[i];
			        		}
			        		
			        		double af;
			        		if(adj_not_found) {
			        			af = Double.valueOf(ac[i])/Double.valueOf(an);
			        		} else {
			        			af = Double.valueOf(ac_adj[i])/Double.valueOf(an_adj);
			        		}
			        		getLoFVariant(chr, pos, ref, id, alt_alleles[i], af,i+1,infos, genotypes);
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
		
		String chr = ""; 
		
		for(LoFVariant lof: lof_statistic)  {
			
			if(!chr.equals(lof.getChr())) {
				for(String key : gene_sample_statistic.keySet()) {
					gene_sample_statistic.get(key).calculate();
				}
				this.extendSampleStatistic();
				this.extendGeneStatistic();
				gene_sample_statistic = new HashMap<>();
				chr = lof.getChr();
				logger.info("Processing chromosome: "+chr);
			}
			
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
					
					if(!(gt.charAt(0)=='0' && gt.charAt(2)=='0')) {
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
		this.extendSampleStatistic();
		this.extendGeneStatistic();
	}
	
	private void extendGeneStatistic() {
		
		for(String gene_sampleId:gene_sample_statistic.keySet()) {
			String gene = gene_sampleId.split("_")[0];
			String sample_id = gene_sampleId.split("_")[1];
			
			GeneInfo gi = gene_statistic.get(gene);

			if (gene_sample_statistic.get(gene_sampleId).isUnaffected()) {
				gi.addUnaffectedSample(sample_id);
			}
			if (sample_statistic.get(sample_id).getPart_LOF_genes().contains(gene)) {
				gi.addAffectedSample(sample_id);
			} 
			
			if (sample_statistic.get(sample_id).getComplete_LOF_genes().contains(gene)) {
				gi.addKOSample(sample_id);
			}
		}
	}
	
	/**
	 * calculate number of full and partial LoFs per sample
	 * depends only on lof_statistic
	 */
	private void generateSampleStatistic() {

		boolean isFull = false;
		
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
	}
	
	/**
	 * extends sample summary by affected and complete inactivated genes
	 * depends on gene_sample_statistic
	 */
	private void extendSampleStatistic() {
		//get complete LOF genes
		for(Entry<String,GeneSampleInfo> e: gene_sample_statistic.entrySet()) {
			int all_trans = gene2transcripts.get(e.getKey().split("_")[0]).size();
			int lof_hom = e.getValue().getHomTrans();
			int lof_het = e.getValue().getHetTrans();
			
			SampleInfo si = sample_statistic.get(e.getKey().split("_")[1]);
			if(all_trans == lof_hom) {//all transcripts are affected by LoF variants
				si.addCompleteLoFGene(e.getKey().split("_")[0]);
				si.addPartLoFGene(e.getKey().split("_")[0]);
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
		
		//average genes per sample
		double avg_aff = 0.0;
		double avg_in = 0.0;
				
		for(String s: sample_statistic.keySet()) {
			SampleInfo si = sample_statistic.get(s);
			avg_aff += si.getPart_LOF_genes().size();
			avg_in += si.getComplete_LOF_genes().size();
		}
				
		avg_aff = avg_aff / sample_statistic.size();
		avg_in = avg_in / sample_statistic.size();
		
		BufferedWriter bw = new BufferedWriter(new FileWriter(outfile));
		bw.write("gene_id\tgene_symbol\tfull\tpartial\tP_LoF_aff\tunaffected_samples\taffected_samples\tp-value_aff\tP_LoF_ko\tko_samples\tp-value_in");

		bw.newLine();
		for(String gene: gene_statistic.keySet()) {

			GeneInfo gi = gene_statistic.get(gene);
			bw.write(gene+"\t"+gi.getSymbol());
			bw.write("\t"+gi.getFullLoFs());
			bw.write("\t"+gi.getPartLoFs());			
			
			double p_lof_aff;
			
			if(gene2background.containsKey(gene)) {
				//deviation factor * overall LoF probability
				p_lof_aff = gene2background.get(gene)[0] * (avg_aff/gene2transcripts.size());
			} else {
				p_lof_aff = 1-gi.getProbUnaffected();
			}
			bw.write("\t"+p_lof_aff);
			
			int unaffected = gi.getUnaffectedSamples().size();
			bw.write("\t"+unaffected);
			int affected = gi.getAffectedSamples().size();
			bw.write("\t"+affected);
			
			int n = affected+unaffected;
			
			double p_val = 1.0;
			if(p_lof_aff>1.0) {
				System.err.println(gene+" "+avg_aff+" "+p_lof_aff);
			} else {
				p_val = 1 - new BinomialDistribution(n, p_lof_aff).cumulativeProbability(affected-1);
			}
			bw.write("\t"+p_val);
			
			double p_lof_in = 0.0;
			if(gene2background.containsKey(gene)) {
				//deviation factor * overall LoF probability
				p_lof_in = gene2background.get(gene)[1] * (avg_in/gene2transcripts.size());
			}
			bw.write("\t"+p_lof_in);
			
			int ko = gi.getKOSamples().size();
			bw.write("\t"+ko);
			
			double p_val_in = 1.0;
			if(p_lof_in >= 1.0) {
				System.err.println(gene+" "+avg_aff+" "+p_lof_in);
			} else {
				p_val_in = 1 - new BinomialDistribution(n, p_lof_in).cumulativeProbability(gi.getKOSamples().size()-1);
			}
			bw.write("\t"+p_val_in);
			
//			double frequency = ((double) affected)/((double) (unaffected+affected));
//			bw.write("\t"+frequency);
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
			int affected = stat.getPart_LOF_genes().size();
			bw.write(s+"\t"+stat.getFullLOFs()+"\t"+stat.getPartLOFs()+"\t"+affected+"\t"+completes.size()+"\t");
			
			if(completes.size() >= 1) {
				String gene = completes.get(0);
				ko_genes.add(gene);
				bw.write(gene +":"+gene_statistic.get(gene).getSymbol());
			}
			
			for(int i = 1; i< completes.size(); i++) {
				String u = completes.get(i);
				ko_genes.add(u);
				bw.write(","+u+":"+gene_statistic.get(u).getSymbol());
			}
			
			String affectedGeneList = "";
			
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
			ko_genes = stat.getComplete_LOF_genes();
			boolean has_parents = false;
			
			if(my_sample_ids.contains(stat.getMatId())) {
				has_parents = true;
				SampleInfo mat = sample_statistic.get(stat.getMatId());
				lof_genes.removeAll(mat.part_LOF_genes);
				ko_genes.removeAll(mat.complete_LOF_genes);
			}
			
			if(my_sample_ids.contains(stat.getPatId())) {
				has_parents = true;
				SampleInfo pat = sample_statistic.get(stat.getPatId());
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
	
	private void writeGeneticBackground(String outfile) throws IOException{
		
		//average LoF/KO genes per sample
		double avg_aff = 0.0;
		double avg_in = 0.0;
		
		for(String s: sample_statistic.keySet()) {
			SampleInfo si = sample_statistic.get(s);
			avg_aff += si.getPart_LOF_genes().size();
			avg_in += si.getComplete_LOF_genes().size();
		}
		
		avg_aff = avg_aff / sample_statistic.size();
		avg_in = avg_in / sample_statistic.size();
		logger.info("Average number of affected genes per sample: " +avg_aff + "/" + avg_in + "(ko genes)");
		
		//all genes / average number of genes per sample
		double norm_aff = (double)gene2transcripts.size()/avg_aff;
		double norm_in = (double)gene2transcripts.size()/avg_in;
		
		if(geneback_header.contains("#trained on:\t"+vcf_file)) {
			logger.info("Genetic background file has yet been trained on: "+vcf_file);
			return;
		}
		// iterate over genes found in analyzed study and integrate numbers in
		// read gene2background
		for (String gene : gene_statistic.keySet()) {
			
			GeneInfo gi = gene_statistic.get(gene);
			Double[] counts = gene2background.get(gene);
			if (counts == null) {
				counts = new Double[3];
				for (int i = 0; i < 3; i++) {
					counts[i] = 0.0;
				}
			}			
			double called_samples = (double)(gi.getAffectedSamples().size() + gi.getUnaffectedSamples().size());
			
			//counts[0] = deviation in background file + deviation in current study weighted by sample number
			counts[0] = ((counts[0] * counts[2]) + (gi.getAffectedSamples().size() * norm_aff))/(called_samples+counts[2]);
			counts[1] = ((counts[1] * counts[2]) + (gi.getKOSamples().size() * norm_in))/(called_samples+counts[2]);
			counts[2] += called_samples;
			gene2background.put(gene, counts);
		}

		
		BufferedWriter bw = new BufferedWriter(new FileWriter(outfile));
		

		for (int i = 0; i < geneback_header.size()-1; i++) {
			bw.write(geneback_header.get(i));
			bw.newLine();
		}
		
		bw.write("#trained on:\t"+vcf_file);
		bw.newLine();
		bw.write("#gene_id\tdeviation_affected\tdeviation_ko\tcalled_samples");
		bw.newLine();
		
		for(String gene:gene2background.keySet()) {
			Double[] counts = gene2background.get(gene);
			bw.write(gene+"\t"+counts[0]+"\t"+counts[1]+"\t"+counts[2]);
			
			bw.newLine();
		}
		bw.close();
	}
	
	private void writeGeneList(HashSet<String> genes, String outfile) throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(outfile));
		for(String gene: genes) {
			
			GeneInfo gi = gene_statistic.get(gene);
			
			bw.write(gene+"\t"+gi.getSymbol());
			
			String samples = "";
			for(String sample_id: gi.getAffectedSamples()) {
				samples += ","+sample_id;
			}
			
			bw.write("\t"+samples.replaceFirst(",", ""));
			bw.newLine();
		}
		bw.close();
	}
}