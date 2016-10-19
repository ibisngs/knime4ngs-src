/**
 *  Copyright (C) 2016 the Knime4NGS contributors.
 *  Website: http://ibisngs.github.io/knime4ngs
 *  
 *  This file is part of the KNIME4NGS KNIME extension.
 *  
 *  The KNIME4NGS extension is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
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
