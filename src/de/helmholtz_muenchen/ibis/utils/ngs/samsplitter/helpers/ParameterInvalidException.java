/**
 * 
 */
package de.helmholtz_muenchen.ibis.utils.ngs.samsplitter.helpers;

/**
 * @author jonathan.hoser
 *
 */
public class ParameterInvalidException extends Exception {

	/**
	 * Generated SVUID
	 */
	private static final long serialVersionUID = -8268895637564821403L;

	/**
	 * Empty proxy towards Exception.
	 * @param message
	 */
	public ParameterInvalidException (String message) 
		{
			super(message);
		}
}
