package de.helmholtz_muenchen.ibis.utils.ngs.frost;

import java.util.ArrayList;

/**
 * @author tanzeem.haque
 * Object for InputScanner
 *
 */
public class InputData {

	private String id;
//	private ArrayList<Position_Type> position_arrList;
	private ArrayList<Integer> positions;
	
	public InputData(String id, ArrayList<Integer> positions) {
		// TODO Auto-generated constructor stub
		setId(id);
		setPositions(positions);
	}
	
//	public InputData(String id, ArrayList<Position_Type> position) {
//		// TODO Auto-generated constructor stub
//		setId(id);
//		setPosition_arrList(position);
//	}

	/**
	 * @return the id
	 */
	protected String getId() {
		return id;
	}

	/**
	 * @param id the id to set
	 */
	protected void setId(String id) {
		this.id = id;
	}

	/**
	 * @return the position
	 */
//	protected ArrayList<Position_Type> getPosition_arrList() {
//		return position_arrList;
//	}
	/**
	 * @param position the position to set
	 */
//	protected void setPosition_arrList(ArrayList<Position_Type> position) {
//		this.position_arrList = position;
//	}

	/**
	 * @return the positions
	 */
	protected ArrayList<Integer> getPositions() {
		return positions;
	}

	/**
	 * @param positions the positions to set
	 */
	protected void setPositions(ArrayList<Integer> positions) {
		this.positions = positions;
	}

	
	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		// TODO Auto-generated method stub
		String data = "";
		for (int i = 0; i < getPositions().size(); i++) {
			data += getId() + "\t" + getPositions().get(i) + "\n";
		}
		return data;
	}

}
