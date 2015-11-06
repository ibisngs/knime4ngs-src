package de.helmholtz_muenchen.ibis.utils.ngs;

import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;

import de.helmholtz_muenchen.ibis.ngs.lofsummary.Annotation;

public interface AnnotationParser {
	
	public String getAnnId();
	public LinkedList<Annotation> getAnnotations(int allele_id, String anno);
	public HashSet<String> getAnnotatedAlleleIds(String anno);
	public HashMap<String, HashSet<Integer>> getGene2AlleleIds (String anno, boolean use_id);
	public HashMap<String, HashSet<String>> getGene2TranscriptIds (int allele_id, String anno);
}