package de.helmholtz_muenchen.ibis.utils.ngs;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Scanner;

public class VCFFile implements Iterator<VCFVariant> {
	
	private String [] vcf_header;
	private HashMap<String, Integer> field_index;
	final String CHR = "CHR";
	final String POS = "POS";
	final String ID = "ID";
	final String REF = "REF";
	final String ALT = "ALT";
	final String QUAL = "QUAL";
	final String FILTER = "FILTER";
	final String INFO = "INFO";
	final String FORMAT = "FORMAT";
	
	private HashMap<String, String> info_header;
	final String INFO_TAG = "##INFO=";
	private Scanner sc;

	
	public VCFFile (String path) throws FileNotFoundException {
		FileInputStream inputStream = new FileInputStream(path);
		this.sc = new Scanner(inputStream, "UTF-8");
		String line = "";
		boolean found_header = false;
		field_index = new HashMap<>();
		info_header = new HashMap<>();
		
		while (!found_header && sc.hasNextLine()) {
			line = sc.nextLine();
			if(line.startsWith("#CHR")) {
				line = line.replaceFirst("#", "");
				vcf_header = line.split("\t");
				for(int i = 0; i < vcf_header.length; i++) {
					String tmp =vcf_header[i];
					if (tmp.contains(CHR)) {
						field_index.put(CHR, i);
					} else if (tmp.equals(POS)) {
						field_index.put(POS, i);
					} else if (tmp.equals(ID)) {
						field_index.put(ID, i);
					} else if (tmp.equals(REF)) {
						field_index.put(REF, i);
					} else if (tmp.equals(ALT)) {
						field_index.put(ALT, i);
					} else if (tmp.equals(QUAL)) {
						field_index.put(QUAL, i);
					} else if (tmp.equals(FILTER)) {
						field_index.put(FILTER, i);
					} else if (tmp.equals(INFO)) {
						field_index.put(INFO, i);
					} else if (tmp.equals(FORMAT)) {
						field_index.put(FORMAT, i);
					}
				}
				
				found_header = true;
			} else if (line.startsWith(INFO_TAG)){
				String field = line.split(INFO_TAG)[1];
				field = field.replaceAll(">", "");
				field = field.replaceAll("<", "");
				String id = field.split(",")[0].split("=")[1];
				info_header.put(id,field);
			}
		}
	}

	@Override
	public boolean hasNext() {
		boolean result = sc.hasNextLine();
		if(!result) {
			sc.close();
		}
		return result;
	}

	@Override
	public VCFVariant next() {
		VCFVariant var = new VCFVariant();
		
		String line = sc.nextLine();
		String fields [] = line.split("\t");
		
		var.setChrom(fields[field_index.get(CHR)]);
		var.setPos(fields[field_index.get(POS)]);
		var.setId(fields[field_index.get(ID)]);
		var.setRef(fields[field_index.get(REF)]);
		var.setAlt(fields[field_index.get(ALT)]);
		var.setQual(fields[field_index.get(QUAL)]);
		var.setFilter(fields[field_index.get(FILTER)]);
		var.setInfo(fields[field_index.get(INFO)]);
		
		if(field_index.containsKey(FORMAT)) {
			var.setFormat(fields[field_index.get(FORMAT)]);
			
			for(int i = 9; i < vcf_header.length; i++) {
				var.addGenotype(vcf_header[i], fields[i]);
			}
		}
		return var;
	}

	@Override
	public void remove() {
		sc.remove();
	}
	
	public String getInfoHeader(String info_id) {
		return this.info_header.get(info_id);
	}
	
	public ArrayList<String> getSampleIds() {
		ArrayList<String> result = new ArrayList<>();
		if(field_index.containsKey(FORMAT)) {
			for(int i = 9; i < vcf_header.length; i++) {
				result.add(vcf_header[i]);
			}
		} else {
			for(int i = 8; i < vcf_header.length; i++) {
				result.add(vcf_header[i]);
			}
		}
		return result;
	}

}
