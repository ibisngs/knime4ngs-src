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
package de.helmholtz_muenchen.ibis.utils.lofs;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Paths;

public class FileInputReader {
	
	private BufferedReader reader;
	
	/**
	 * creates a buffered reader
	 * @param path: file to read
	 */
	public FileInputReader(String path) throws IOException{
		Charset c = Charset.forName("UTF-8");
		reader = Files.newBufferedReader(Paths.get(path),c);

	}
	
	/**
	 * reads one line from the file
	 * @return
	 */
	public String read(){

			String res;
			try {
				res = reader.readLine();
				return res;
			} catch (IOException e) {
				e.printStackTrace();
			}

		return "";
	}
	
	/**
	 * closes the writer
	 * @throws IOException
	 */
	public void closer() throws IOException{
			reader.close();
	}


}
