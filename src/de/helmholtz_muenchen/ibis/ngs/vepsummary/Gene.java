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

package de.helmholtz_muenchen.ibis.ngs.vepsummary;

import java.util.HashSet;

/**
* 
* @author Tim Jeske
*  
*/

public class Gene {
	
	private String id;
	private String symbol;
	private HashSet<String> transcript_ids;
	
	public Gene(String id) {
		this.id = id;
		this.symbol = "NA";
		this.transcript_ids = new HashSet<>();
	}
	
	public Gene(String id, String symbol) {
		this.id = id;
		this.symbol = symbol;
		this.transcript_ids = new HashSet<>();
	}
	
	public void addTranscript(String t) {
		this.transcript_ids.add(t);
	}
	
	public HashSet<String> getTranscripts() {
		return this.transcript_ids;
	}
	
	public String getSymbol()  {
		return this.symbol;
	}
	
	public String getId() {
		return this.id;
	}
	
	public String toString() {
		String res = id+"\t"+symbol+"\t";
		String transcripts = "";
		for(String t: transcript_ids) {
			transcripts += ","+t;
		}
		transcripts = transcripts.replaceFirst(",", "");
		
		return res+transcripts;
	}

}
