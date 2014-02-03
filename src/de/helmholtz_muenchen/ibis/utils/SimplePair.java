package de.helmholtz_muenchen.ibis.utils;

import org.apache.commons.lang3.tuple.Pair;


public class SimplePair<L,R> extends Pair<L,R> {
	public static final String SEP = "=";
	private static final long serialVersionUID = 7909959941935253757L;
	private L left;
	private R right;
	
	public SimplePair(L l, R r){
		left  = l;
		right = r;
	}
	
	public SimplePair(){
		left  = null;
		right = null;
	}

	@Override
	public L getLeft() {
		return(left);
	}

	@Override
	public R getRight() {
		return(right);
	}

	@Override
	public R setValue(R value) {
		R oldR = right;
		right = value;
		return(oldR);
	}

	@Override
	public String toString(){
		return(left + SEP + right);
	}


}
