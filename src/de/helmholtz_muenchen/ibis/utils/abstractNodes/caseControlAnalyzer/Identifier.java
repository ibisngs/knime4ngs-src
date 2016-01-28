package de.helmholtz_muenchen.ibis.utils.abstractNodes.caseControlAnalyzer;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

import de.helmholtz_muenchen.ibis.utils.ngs.BioEntity;

public interface Identifier {
	
	List<String> getMappings(String identifier);
	
	public class GeneIdentifier implements Identifier {

		@Override
		public List<String> getMappings(String identifier) {
			List<String> res = new ArrayList<>();
			res.add(identifier.split("_",2)[1]); //returns ensemblid_genesymbol
			return res;
		}
	}
	
	public class EntityIdentifier implements Identifier {

		BioEntity e;
		
		public EntityIdentifier(BioEntity e) {
			this.e = e;
		}
		
		@Override
		public List<String> getMappings(String identifier) {
			List<String> res = new ArrayList<String>();
			
			switch(e) {
			case GENE_ID:
				res.add(identifier.split("_")[1].toUpperCase()); //returns ensemblid
				break;

			case GENE_SYMBOL:
				res.add(identifier.split("_")[2].toUpperCase());
				break;
				
			case TRANSCRIPT_ID:	
				res.add(identifier.split("_")[0].toUpperCase()); //returns transcriptid
				break;
			}
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
