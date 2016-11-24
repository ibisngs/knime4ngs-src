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
package de.helmholtz_muenchen.ibis.utils;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.util.Iterator;
import java.util.NoSuchElementException;

public class MMapFile {

	public class MMapIterator implements Iterator<String> {
		private int offset;

		public MMapIterator(int offset) {
			this.offset = offset;
		}

		public boolean hasNext() {
			return offset < cb.limit();
		}

		public String next() {
			ByteArrayOutputStream sb = new ByteArrayOutputStream();
			if (offset >= cb.limit())
				throw new NoSuchElementException();
			for (; offset < cb.limit(); offset++) {
				byte c = (cb.get(offset));
				if (c == '\n') {
					offset++;
					break;
				}
				if (c != '\r') {
					sb.write(c);
				}

			}
			try {
				return sb.toString("UTF-8");
			} catch (UnsupportedEncodingException e) {
			}
			return sb.toString();
		}

		public void remove() {

		}
	}

	private ByteBuffer cb;
	long size;
	private long numLines = -1;

	public MMapFile(String file) throws FileNotFoundException, IOException {
		FileInputStream fileInputStream = new FileInputStream(new File(file));
		FileChannel fc = fileInputStream.getChannel();
		size = fc.size();
		cb = fc.map(FileChannel.MapMode.READ_ONLY, 0, fc.size());
		fileInputStream.close();
	}

	public long getNumLines() {
		if (numLines != -1)
			return numLines; // cache number of lines
		long cnt = 0;
		for (int i = 0; i < size; i++) {
			if (cb.get(i) == '\n')
				cnt++;
		}
		numLines = cnt;
		return cnt;
	}

	public Iterator<String> tail(long lines) {
		long cnt = 0;
		long i = 0;
		for (i = size - 1; i >= 0; i--) {
			if (cb.get((int) i) == '\n') {
				cnt++;
				if (cnt == lines + 1)
					break;
			}
		}
		return new MMapIterator((int) i + 1);
	}

	public Iterator<String> head() {
		return new MMapIterator(0);
	}

	/*static public void main(String[] args) {
		try {
			Iterator<String> it = new MMapFile("/test.txt").head();
			while (it.hasNext()) {
				System.out.println(it.next());
			}
		} catch (Exception e) {

		}

		System.out.println();

		try {
			Iterator<String> it = new MMapFile("/test.txt").tail(2);
			while (it.hasNext()) {
				System.out.println(it.next());
			}
		} catch (Exception e) {

		}

		System.out.println();

		try {
			System.out.println("lines: "
					+ new MMapFile("/test.txt").getNumLines());
		} catch (Exception e) {

		}

	}*/

}
