package de.helmholtz_muenchen.ibis.utils.ngs.frost;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;

/**
 * @author tanzeem.haque
 * Class to write all info outputs
 *
 */

public class OutputWriters {

	private FastaReader fs;
	private InputScanner in;
	private Trio_Simulator trio;
	
	public OutputWriters(FastaReader fs, InputScanner in, Trio_Simulator trio) {
		// TODO Auto-generated constructor stub
		setFs(fs);
		setIn(in);
		setTrio(trio);
	}
	
	protected void write_simple_string(String fileName, String s) {
		File file = new File(fileName);
		try {
			BufferedWriter bw = new BufferedWriter( new FileWriter(file));
			bw.write(s + "\n");	
			bw.close();
		} catch (Exception e) {
			System.err.println("Error: " + e.getMessage());
			e.printStackTrace();
		}		
	}
	
	protected void write_InputData_arrList(String fileName, ArrayList<InputData> iData) {
		File file = new File(fileName);
		try {
			BufferedWriter bw = new BufferedWriter( new FileWriter(file));
			for (int i = 0; i < iData.size(); i++) {
				for (int j = 0; j < iData.get(i).getPositions().size(); j++) {
					bw.write(iData.get(i).getId() + "\t" + iData.get(i).getPositions().get(j) + "\n");
				}
			}	
			bw.close();
		} catch (Exception e) {
			System.err.println("Error: " + e.getMessage());
			e.printStackTrace();
		}		
	}
	
	protected void write_vcf(String fileName, ArrayList<String> ID, 
			ArrayList<Integer> POS, ArrayList<String> prettyMutation) {
		// TODO Auto-generated method stub
		File file = new File(fileName);
		try {
			BufferedWriter bw = new BufferedWriter( new FileWriter(file));
			for (int i = 0; i < prettyMutation.size(); i++) {
				bw.write(ID.get(i) + "\t" + POS.get(i) + "\t" + prettyMutation.get(i) + "\n");
			}	
			bw.close();
		} catch (Exception e) {
			System.err.println("Error: " + e.getMessage());
			e.printStackTrace();
		}		
	}

	/**
	 * @return the fs
	 */
	public FastaReader getFs() {
		return fs;
	}

	/**
	 * @param fs the fs to set
	 */
	public void setFs(FastaReader fs) {
		this.fs = fs;
	}

	/**
	 * @return the in
	 */
	public InputScanner getIn() {
		return in;
	}

	/**
	 * @param in the in to set
	 */
	public void setIn(InputScanner in) {
		this.in = in;
	}

	/**
	 * @return the family
	 */
	public Trio_Simulator getTrio() {
		return trio;
	}

	/**
	 * @param trio the family to set
	 */
	public void setTrio(Trio_Simulator trio) {
		this.trio = trio;
	}

}
