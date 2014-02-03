package de.helmholtz_muenchen.ibis.metabolomics.graphicalModels;

import java.util.ArrayList;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;

public class ParamUtils {
//	/**
//	 * Parse String to HashMap
//	 * @param s the String which shall be parsed to the hashMap
//	 * @param sep1 separator between pairs
//	 * @param sep2 separator between keys and values
//	 * @return the parsed HashMap
//	 * @throws Exception if a pair contains != 2 entries
//	 */
//	public static HashMap<String,String> HashMapFromString(String s, String sep1, String sep2) throws Exception{
//		HashMap<String,String> result = new HashMap<String,String>();
//		if( s == null || s.equals("")){
//			return(result);
//		}
//		String[] pairs = s.split(sep1);
//
//		for(int p=0; p<pairs.length; p++){
//			if( pairs[p] != null && pairs[p] != "" ){
//				String[] pair = pairs[p].split(sep2);
//				if(pair.length != 2){
//					throw(new Exception("Can't parse String '" + pairs[p] + "' to Pair (separated by '" + sep2 + "')"));
//				}
//				result.put(pair[0], pair[1]);
//			}
//		}
//
//		return(result);
//	}
//	public static Map<String, Map<String,String>> HashMapFromString_2dim(String s, String sep1, String sep2, String sep3, String sep4) throws Exception{
//		Map<String, Map<String,String>> result = new HashMap<String,Map<String,String>>();
//
//		String[] pairs = s.split(sep1);
//
//		for(int p=0; p<pairs.length; p++){
//			if( pairs[p] != null && pairs[p] != "" ){
//				String[] pair = pairs[p].split(sep2);
//				if(pair.length != 2){
//					throw(new Exception("Can't parse String '" + pairs[p] + "' to Pair (separated by '" + sep2 + "')"));
//				}
//				result.put(pair[0], HashMapFromString(pair[1], sep3, sep4));
//			}
//		}
//
//		return(result);
//	}
//	public static Map<String, Map<String,String>> HashMapFromString_2dim(String s) throws Exception{
//		return(HashMapFromString_2dim(s, ";", ":", ",", "="));
//	}
//
//
//
//
//
//
//	/**
//	 * Parse HashMap to String
//	 * @param h the HashMap which shall be parsed to the String
//	 * @param sep1 separator between pairs
//	 * @param sep2 separator between keys and values
//	 * @return the parsed String
//	 */
//	public static String HashMapToString(HashMap<String, String> h, String sep1, String sep2){
//		String result = "";
//		Iterator<String> it = h.keySet().iterator();
//		while(it.hasNext()){
//			String k = it.next();
//			result += k + sep2 + h.get(k) + sep1;
//		}
//		result.replace(sep1+"$", "");
//		return(result);
//	}
//	public static String HashMapToString_2dim(HashMap<String, HashMap<String,String>> h, String sep1, String sep2, String sep3, String sep4){
//		String result = "";
//		Iterator<String> it = h.keySet().iterator();
//		while(it.hasNext()){
//			String k = it.next();
//			result += k + sep2 + HashMapToString(h.get(k),sep3,sep4) + sep1;
//		}
//		result.replace(sep1+"$", "");
//		return(result);
//	}
//	public static String HashMapToString_2dim(HashMap<String, HashMap<String,String>> h){
//		return(HashMapToString_2dim(h, ";", ":", ",", "="));
//	}
//
//
//
//
//
//
//	public static <T> HashMap<T,T> arraysToHashMap(T[] keys, T[] values) throws Exception{
//		HashMap<T,T> map = new HashMap<T,T>();
//		if(keys.length != values.length){
//			throw(new Exception("Can't parse arrays of different length to HashMap"));
//		}
//		for(int index = 0; index < keys.length; index++){
//			map.put(keys[index], values[index]);
//		}
//
//		return(map);
//	}


	static String paramsToString(ArrayList<Pair<String,String>> params, String paramsSep, String keyValueSep){
		String result = "";
		if(params == null || params.size() == 0){
			return(result);
		}
		for(int i=0; i<params.size(); i++){
			result += params.get(i).getKey() + keyValueSep + params.get(i).getValue() + paramsSep;
		}
		return(result.substring(0, result.length() - paramsSep.length()));
	}
	
	
	
	static ArrayList<Pair<String,String>> stringToParams( String params, String paramsSep, String keyValueSep) throws Exception{
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
