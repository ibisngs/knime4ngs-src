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

import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;

import de.helmholtz_muenchen.ibis.ngs.vepsummary.Annotation;

public interface AnnotationParser {
	
	public String getAnnId();
	public double getExACAF(int allele_id, String anno);
	public LinkedList<Annotation> getAnnotations(int allele_id, String anno);
	public HashSet<String> getAnnotatedAlleleIds(String anno);
	public HashMap<String, HashSet<Integer>> getEntity2AlleleIds (String anno, BioEntity entitiy);
	public HashMap<String, HashSet<String>> getGene2TranscriptIds (int allele_id, String anno);
}

