package de.helmholtz_muenchen.ibis.utils.ngs.frost;


public class Strand {
	
	private int start;
	private int end;
	private String stream;
	
	public Strand(int start, int end, String stream) {
		setStart(start);
		setEnd(end);
		setStream(stream);
	}

	

	public String getStream() {
		return stream;
	}

	public void setStream(String stream) {
		this.stream = stream;
	}



	public int getStart() {
		return start;
	}



	public void setStart(int start) {
		this.start = start;
	}



	public int getEnd() {
		return end;
	}



	public void setEnd(int end) {
		this.end = end;
	}
	
}
