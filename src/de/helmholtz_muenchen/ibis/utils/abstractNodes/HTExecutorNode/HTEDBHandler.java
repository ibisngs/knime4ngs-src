package de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode;

import java.sql.*;
import java.util.HashSet;

public class HTEDBHandler {
	
	private Connection con;
	
	private final String HTEXECUTION = "CREATE TABLE HTExecution " +
			"(exec_id INTEGER PRIMARY KEY AUTOINCREMENT, " +
			"klock_string BLOB, " +
			"node_name VARCHAR(255), " +
			"host_name VARCHAR(255), " +
			"node_time DATETIME DEFAULT CURRENT_TIMESTAMP, " +
			"node_successful_finished int DEFAULT 0, " +
			"node_threshold int , "+
			"node_count int DEFAULT 0)";
	
	private final String HTERROR = "CREATE TABLE HTError_History " +
			"(err_id INTEGER PRIMARY KEY AUTOINCREMENT, " +
			"exec_id int references HTExecution(exec_id), " +
			"time DATETIME DEFAULT CURRENT_TIMESTAMP, " +
			"err_msg BLOB)";
	
	public HTEDBHandler(String file) throws SQLException {
		try {
			Class.forName("org.sqlite.JDBC").newInstance();
			this.con = DriverManager.getConnection("jdbc:sqlite:"+file);
		} catch (InstantiationException | IllegalAccessException
				| ClassNotFoundException e) {
			e.printStackTrace();
			System.err.println("Failure while loading jdbc driver: "+e.getMessage());
		}
	}
	
	public boolean checkSchema() throws SQLException {
		
		HashSet<String> tables = new HashSet<>();
		tables.add("HTExecution");
		tables.add("HTError_History");
		Statement stmt = con.createStatement();
		String query = "SELECT name,sql FROM sqlite_master WHERE type = 'table'";
		ResultSet rs = stmt.executeQuery(query);
		while(rs.next()) {
			String name = rs.getString("name");
			String def = rs.getString("sql");
			if(def.equals(HTERROR) || def.equals(HTEXECUTION)) {
				tables.remove(name);
			}
		}
		
		return tables.isEmpty();
	}
	
	public void createDB() throws SQLException {
		Statement stmt = con.createStatement();
		stmt.executeUpdate(HTEXECUTION);

		stmt.executeUpdate(HTERROR);
		stmt.close();
	}
	
	public int insertNewHTExecution(String command, String name, String host, int node_threshold) {
		String query = "INSERT INTO HTExecution (klock_string, node_name, host_name, node_threshold) VALUES (?,?,?,?)";
		PreparedStatement preStmt;
		ResultSet rs;
		int id = -1;
		try {
			preStmt = con.prepareStatement(query,Statement.RETURN_GENERATED_KEYS);
			preStmt.setString(1, command);
			preStmt.setString(2, name);
			preStmt.setString(3, host);
			preStmt.setInt(4, node_threshold);
			preStmt.execute();
			rs = preStmt.getGeneratedKeys();
			if(rs.next()) {
				id = rs.getInt(1);
			}
		} catch (SQLException e) {
			System.err.println("New HTExecution could not be inserted: "+e.getMessage());
		}
		return id;
	}
	
	public void updateCount(int exec_id, int count) {
		String query = "UPDATE HTExecution SET node_count = ? WHERE exec_id = ?";
		PreparedStatement preStmt;
		try {
			preStmt = con.prepareStatement(query);
			preStmt.setInt(1, count);
			preStmt.setInt(2, exec_id);
			preStmt.execute();
		} catch (SQLException e) {
			System.err.println("Count could not be updated for node execution with id "+exec_id+": "+e.getMessage());
		}
	}
	
	public void writeSuccess(int exec_id) {
		String query = "UPDATE HTExecution SET node_successful_finished = 1 WHERE exec_id = ?";
		PreparedStatement preStmt;
		try {
			preStmt = con.prepareStatement(query);
			preStmt.setInt(1, exec_id);
			preStmt.execute();
		} catch (SQLException e) {
			System.err.println("Successful could not be updated for node execution with id "+exec_id+": "+e.getMessage());
		}
	}
	
	public void writeError(int exec_id, String errmsg) {
		String query = "INSERT INTO HTError_History (exec_id, err_msg) VALUES (?,?);";
		PreparedStatement preStmt;
		try {
			preStmt = con.prepareStatement(query);
			preStmt.setInt(1, exec_id);
			preStmt.setString(2, errmsg);
			preStmt.execute();
		} catch (SQLException e) {
			System.err.println("Error message could not be inserted for node execution with id "+exec_id+": "+e.getMessage());
		}
	}
	
	public void closeConnection() {
		try {
			this.con.close();
		} catch (SQLException e) {
			System.err.println("Connection to database couldn't be closed correctly: "+e.getMessage());
		}
	}
}
