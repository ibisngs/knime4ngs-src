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

/**
* 
*  * @author Tim Jeske
*  
*/

public class Annotation {
	
	private String gene_id;
	private String transcript_id;
	private String consequence;
	
	public Annotation(String gene_id, String transcript_id, String consequence) {
		this.gene_id = gene_id;
		this.transcript_id = transcript_id;
		this.consequence = consequence;
	}
	
	public String toString() {
		return "("+gene_id+"|"+transcript_id+"|"+consequence+")";
	}

}
