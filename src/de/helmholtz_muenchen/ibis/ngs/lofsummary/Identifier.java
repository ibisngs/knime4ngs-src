package de.helmholtz_muenchen.ibis.ngs.lofsummary;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

public interface Identifier {
	
	List<String> getMappings(String identifier);
	
	public class GeneIdentifier implements Identifier {

		@Override
		public List<String> getMappings(String identifier) {
			List<String> res = new ArrayList<>();
			res.add(identifier.split("_",2)[1]); //id_symbol
			return res;
		}
	}
	
	public class GeneSymbolIdentifier implements Identifier {

		@Override
		public List<String> getMappings(String identifier) {
			List<String> res = new ArrayList<String>();
			res.add(identifier.split("_")[2].toUpperCase());
			return res;
		}
	}
	
	public class GeneSetIdentifier implements Identifier {

		HashMap <String, HashSet<String>> gene2sets;
		
		public GeneSetIdentifier(HashMap <String, HashSet<String>> gene2sets) {
			this.gene2sets = gene2sets;
		}
		
		@Override
		public List<String> getMappings(String identifier) {
			String gene = identifier.split("_")[2].toUpperCase();
			List<String> res = new ArrayList<String>();
			if(gene2sets.containsKey(gene)) {
				res.addAll(gene2sets.get(gene));
			}
			return res;
		}
	}

}
