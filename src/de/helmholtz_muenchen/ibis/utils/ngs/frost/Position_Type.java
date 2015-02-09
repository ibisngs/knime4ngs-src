package de.helmholtz_muenchen.ibis.utils.ngs.frost;

/**
 * @author tanzeem.haque
 * Object for Mutation
 *
 */
public class Position_Type implements Comparable<Position_Type>{

	private int pos;
	private char type;
	
	public Position_Type(int pos, char type) {
		// TODO Auto-generated constructor stub
		setPos(pos);
		setType(type);
	}

	/**
	 * @return the pos
	 */
	public int getPos() {
		return pos;
	}

	/**
	 * @param pos the pos to set
	 */
	public void setPos(int pos) {
		this.pos = pos;
	}

	/**
	 * @return the type
	 */
	public char getType() {
		return type;
	}

	/**
	 * @param type the type to set
	 */
	public void setType(char type) {
		this.type = type;
	}

	@Override
	public int compareTo(Position_Type position_type) {
		// TODO Auto-generated method stub
		int comparePosition = ((Position_Type) position_type).getPos(); 
		return this.pos - comparePosition;
	}

}
