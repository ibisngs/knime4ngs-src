/**
 *  Copyright (C) 2016 the Knime4NGS contributors.
 *  Website: http://ibisngs.github.io/knime4ngs
 *  
 *  This file is part of the KNIME4NGS KNIME extension.
 *  
 *  The KNIME4NGS extension is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package de.helmholtz_muenchen.ibis.utils.ngs;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Scanner;

import org.apache.xmlbeans.SystemProperties;

public class VCFFile implements Iterator<VCFVariant> {
	
	private String [] chr_header;
	private String header;
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
		StringBuilder header_lines = new StringBuilder();
		field_index = new HashMap<>();
		info_header = new HashMap<>();
		
		while (!found_header && sc.hasNextLine()) {
			line = sc.nextLine();
			if(line.startsWith("#")) {
				header_lines.append(line.trim()+SystemProperties.getProperty("line.separator"));
			}
			if(line.startsWith("#CHR")) {
				line = line.replaceFirst("#", "");
				chr_header = line.split("\t");
				for(int i = 0; i < chr_header.length; i++) {
					String tmp =chr_header[i];
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
		
		this.header = header_lines.toString();
	}

	public String getCompleteHeader() {
		return this.header;
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
			
			for(int i = 9; i < chr_header.length; i++) {
				var.addGenotype(chr_header[i], fields[i]);
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
			for(int i = 9; i < chr_header.length; i++) {
				result.add(chr_header[i]);
			}
		} else {
			for(int i = 8; i < chr_header.length; i++) {
				result.add(chr_header[i]);
			}
		}
		return result;
	}

}
