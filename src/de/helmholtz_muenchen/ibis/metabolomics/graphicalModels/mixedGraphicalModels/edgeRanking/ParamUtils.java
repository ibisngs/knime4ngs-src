package de.helmholtz_muenchen.ibis.metabolomics.graphicalModels.mixedGraphicalModels.edgeRanking;

import java.util.ArrayList;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;

public class ParamUtils {
	public static String paramsToString(ArrayList<Pair<String,String>> params, String paramsSep, String keyValueSep){
		String result = "";
		if(params == null || params.size() == 0){
			return(result);
		}
		for(int i=0; i<params.size(); i++){
			result += params.get(i).getKey() + keyValueSep + params.get(i).getValue() + paramsSep;
		}
		return(result.substring(0, result.length() - paramsSep.length()));
	}
	
	
	
	public static ArrayList<Pair<String,String>> stringToParams( String params, String paramsSep, String keyValueSep) throws Exception{
		ArrayList<Pair<String,String>> result = new ArrayList<Pair<String,String>>();
		if(params==null || params.matches("^\\s*$")){
			return(result);
		}

		String[] pairs = params.split(paramsSep);
		for(int p=0; p<pairs.length; p++){
			if(pairs[p].matches("^\\s*$")){
				continue;
			}
			String[] pair = pairs[p].split(keyValueSep);
			if(pair.length != 2){
				throw(new Exception("Can't split "+ pairs[p] +" into key"+keyValueSep+" pair"));
			}
			result.add(new ImmutablePair<String,String>(pair[0],pair[1]));
		}
		return(result);
	}
	
	
}
