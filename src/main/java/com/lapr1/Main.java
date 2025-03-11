package com.lapr1;

import java.io.*;
import java.util.*;

public class Main {
	//caminho do gnuplot
	public final static String GNUPLOT = "c:\\Program Files\\gnuplot\\bin\\gnuplot.exe";

	public static void main (String[] args) throws IOException, InterruptedException {

		Parametros parametros;

		if (args.length != 0 && args.length!=9) {
			System.out.println("O número de argumentos é inválido");
			return;
		}

		if (args.length==0){
			//modo interactivo
			parametros = askParameters();
			List<String[]> parametrosModelo = readFile(parametros.ficheiroParametros);
			for (int i = 0; i < parametrosModelo.size(); i++){

				String namePerson = parametrosModelo.get(i)[0];

				if (Objects.equals(namePerson, parametros.namePerson)) {
					parametros.beta = Double.valueOf(parametrosModelo.get(i)[1]);
					parametros.gama = Double.valueOf(parametrosModelo.get(i)[2]);
					parametros.ro = Double.valueOf(parametrosModelo.get(i)[3]);
					parametros.alfa = Double.valueOf(parametrosModelo.get(i)[4]);
				}
			}

			//beta;gama;ro;alfa
			double beta = parametros.beta;
			double gama = parametros.gama;
			double ro   = parametros.ro;
			double alfa = parametros.alfa;

			double N = parametros.populationN;
			double S = N - 1;
			double I = 1;
			double R = 0;

			if (parametros.metodo.equals("1")) {
				double [][] resultsEuler = euler(S, I, R, parametros.stepH, parametros.daysn, beta, alfa, gama, ro);
				writeFile(resultsEuler, parametros.ficheiroResultados);
			}else{
				double[][] resultsRk4 = rungeKutta4(S, I, R, parametros.stepH, parametros.daysn, beta, alfa, gama, ro);
				writeFile(resultsRk4, parametros.ficheiroResultados);
			}

		}else{
			
			parametros = readCommandLine(args);

			List<String[]> list = readFile(parametros.ficheiroParametros /*~"src/main/Files/ficheiroSIR.csv"*/);

			List<Parametros> parametrosList = transformValues(list);

			double N = parametros.populationN;
			double S = N -1;
			double I = 1;
			double R = 0;

			for (Parametros parametrosEntry : parametrosList) {

				double beta = parametrosEntry.beta;
				double gama = parametrosEntry.gama;
				double ro   = parametrosEntry.ro;
				double alfa = parametrosEntry.alfa;
				String namePerson = parametrosEntry.namePerson;
				parametros.ficheiroResultados = namePerson + parametros.ficheiroResultadosNome;

				if (parametros.metodo.equals("1")) {
					double [][] resultsEuler = euler(S, I, R, parametros.stepH, parametros.daysn, beta, alfa, gama, ro);
					writeFile(resultsEuler, parametros.ficheiroResultados);
				}else{
					double[][] resultsRk4 = rungeKutta4(S, I, R, parametros.stepH, parametros.daysn, beta, alfa, gama, ro);
					writeFile(resultsRk4, parametros.ficheiroResultados);
				}

			}

		}
		
		if (args.length == 0){
			System.out.println("Pretende gerar gráfico (S/N) ?");
			Scanner sc = new Scanner(System.in);

			boolean gerarGrafico = sc.nextLine().equalsIgnoreCase("s");

			if (gerarGrafico)
				gerarGrafico(parametros, parametros.ficheiroResultados.substring(0,parametros.ficheiroResultados.length()-4));
			
		}
	}

	private static void gerarGrafico(Parametros parametros, String ficheiro) throws IOException, InterruptedException {
		//primeiro geraar o ficheiro com o script do gnuplot

		FileWriter fw = new FileWriter(ficheiro+".txt");
		fw.write("set title 'Distribuição de Noticias Falsas'\n");
		fw.write("set datafile separator \";\"\n");
		fw.write("set xlabel 'Dias'\n");
		fw.write("set ylabel 'População'\n");
		fw.write("set key autotitle columnhead\n");
		fw.write("set terminal pngcairo\n");
		fw.write("set output '"+ficheiro+".png'\n");
		fw.write("set decimalsign ','\n");
		fw.write("set xrange[0:" + parametros.daysn + "]\n");
		fw.write("set yrange[0:" + (parametros.populationN + 50) + "]\n");
		//fw.write("set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 pi -1 ps 1.5\n");

		fw.write("plot '"+ficheiro+".csv' using 0:2 with lines, '"
				+ficheiro+".csv' using 0:3 with lines, '"
				+ficheiro+".csv' using 0:4 with lines, '"
				+ficheiro+".csv' using 0:5 with lines\n");

		fw.write("set output\n");
		fw.close();

		String[] commands = new String[2];
		commands[0] = GNUPLOT;
		commands[1] = ficheiro+".txt";


		Process p = Runtime.getRuntime().exec(commands);

		p.waitFor();

	}


	public static double fS(double S, double I, double beta) {
		return (-beta) * S * I;
	}

	public static double fI(double S, double I, double R, double beta, double alfa, double gama, double ro) {
		return ro * beta * S * I - gama * I + alfa * R;
	}

	public static double fR(double S, double I, double R, double beta, double alfa, double gama, double ro) {
		return gama * I - alfa * R + (1 - ro) * beta * S * I;
	}

	public static double[][] euler(double S, double I, double R, double h, int n, double beta, double alfa, double gama, double ro) {

		int i = 0;
		double[][] resultsEuler = new double[n][5];
		resultsEuler[i][0] = i;
		resultsEuler[i][1] = S;
		resultsEuler[i][2] = I;
		resultsEuler[i][3] = R;
		resultsEuler[i][4] = S+I+R;

		for (i = 1; i < n; i++) {
			for (float j = 0; j < 1; j+=h) {

				double Sn = S + h * fS(S, I, beta);
				double In = I + h * fI(S, I, R, beta, alfa, gama, ro);
				double Rn = R + h * fR(S, I, R, beta, alfa, gama, ro);
				S = Sn;
				I = In;
				R = Rn;

			}
			resultsEuler[i][0] = i;
			resultsEuler[i][1] = S;
			resultsEuler[i][2] = I;
			resultsEuler[i][3] = R;
			resultsEuler[i][4] = S+I+R;



		}

		return resultsEuler;
	}

	public static double [][] rungeKutta4 (double S, double I, double R, double h, int n, double beta, double alfa, double gama, double ro) {

		int i = 0;
		double [][] resultsRk4 = new double [n][5];
		resultsRk4[i][0] = i;
		resultsRk4[i][1] = S;
		resultsRk4[i][2] = I;
		resultsRk4[i][3] = R;
		resultsRk4[i][4] = S+I+R;

		for (i=1; i < n; i++) {
			for (float j = 0; j < 1; j+=h) {

				double k1S = h * fS(S, I, beta);
				double k1I = h * fI(S, I, R, beta, alfa, gama, ro);
				double k1R = h * fR(S, I, R, beta, alfa, gama, ro);

				double k2S = h * fS(S + k1S/2, I + k1I/2, beta);
				double k2I = h * fI(S + k1S/2, I + k1I/2, R + k1R/2, beta, alfa, gama, ro);
				double k2R = h * fR(S + k1S/2, I + k1I/2, R + k1R/2, beta, alfa, gama, ro);

				double k3S = h * fS(S + k2S/2, I + k2I/2, beta);
				double k3I = h * fI(S + k2S/2, I + k2I/2, R + k2R/2, beta, alfa, gama, ro);
				double k3R = h * fR(S + k2S/2, I + k2I/2, R + k2R/2, beta, alfa, gama, ro);

				double k4S = h * fS(S + k3S, I + k3I, beta);
				double k4I = h * fI(S + k3S, I + k3I, R + k3R, beta, alfa, gama, ro);
				double k4R = h * fR(S + k3S, I + k3I, R + k3R, beta, alfa, gama, ro);

				double kS = (k1S + (2 * k2S) + (2 * k3S) + k4S) / 6;
				double kI = (k1I + (2 * k2I) + (2 * k3I) + k4I) / 6;
				double kR = (k1R + (2 * k2R) + (2 * k3R) + k4R) / 6;

				double Sn = S + kS;
				double In = I + kI;
				Double Rn = R + kR;

				S = Sn;
				I = In;
				R = Rn;


			}

			resultsRk4[i][0] = i;
			resultsRk4[i][1] = S;
			resultsRk4[i][2] = I;
			resultsRk4[i][3] = R;
			resultsRk4[i][4] = S+I+R;

		}
		return resultsRk4;
	}


	public static Parametros askParameters() {
		Parametros parametros = new Parametros();
		Scanner sc = new Scanner(System.in);

		System.out.println("Introduza o Ficheiro de Parâmetros");
		parametros.ficheiroParametros = sc.nextLine();

		System.out.println("Introduza o nome da pessoa para a leitura de parametros");
		parametros.namePerson = sc.nextLine();

		boolean invalido;
		do {
			invalido = false;
			System.out.println("Introduza o Método: 1-Euler ou 2-Runge Kutta");
			parametros.metodo = sc.nextLine();
			if (!parametros.metodo.equals("1") && !parametros.metodo.equals("2") ) {
				invalido = true;
				System.out.println("O Parâmetro -m é inválido, tem de ser 1 ou 2");
			}
		}while(invalido);

		do {
			invalido = false;
			System.out.println("Introduza o valor do Passo de Integracão");
			parametros.stepH = sc.nextDouble();

			if (parametros.stepH <= 0 || parametros.stepH >= 1) {
				invalido = true;
				System.out.println("O Parametro -p é inválido, tem de ser entre 0 e 1");
			}
		}while(invalido);

		System.out.println("Introduza o Tamanho da População");
		parametros.populationN = sc.nextInt();

		System.out.println("Introduza o Número de Dias");
		parametros.daysn = sc.nextInt();
		sc.nextLine(); // para ler a mudanca de linha

		System.out.println("Introduza o Nome do ficheiro de resultados");
		parametros.ficheiroResultados = sc.nextLine();

		return parametros;
	}

	public static Parametros readCommandLine(String[] args) {
		Parametros parametros = new Parametros();
		
		parametros.ficheiroParametros = args[0];
		if (!args[1].equals("-m")){
			System.out.println("Falta o parametro -m");
			return null;
		}
		parametros.metodo=args[2];

		if (!parametros.metodo.equals("1") && !parametros.metodo.equals("2")){
			System.out.println("O parâmetro -m tem de ser 1-Euler ou 2-Range Kutta");
			return null;
		}

		if (!args[3].equals("-p")){
			System.out.println("Falta o parâmetro -p");
			return null;
		}
		try {
			parametros.stepH = Double.parseDouble(args[4]);

			if (parametros.stepH <=0 || parametros.stepH >= 1){
				System.out.println("O parâmetro -p é inválido, tem de ser entre 0 e 1");
				return null;
			}

		} catch(NumberFormatException ex){
			System.out.println("O parâmetro -p é inválido");
			return null;
		}

		if (!args[5].equals("-t")){
			System.out.println("Falta o parâmetro -t");
			return null;
		}
		try {
			parametros.populationN = Integer.parseInt(args[6]);
		} catch(NumberFormatException ex){
			System.out.println("O parâmetro -t é inválido");
			return null;
		}

		if (!args[7].equals("-d")){
			System.out.println("Falta o parâmetro -d");
			return null;
		}
		try {
			parametros.daysn = Integer.parseInt(args[8]);
		}catch(NumberFormatException ex){
			System.out.println("O Parâmetro -d é inválido");
			return null;
		}

		args[4] = args[4].replace(".", "");
		parametros.ficheiroResultadosNome = "m"+args[2]+"p"+args[4]+"t"+args[6]+"d"+args[8]+".csv";
		return parametros;
	}

	public static void writeFile (double[][] matriz, String filename) throws IOException {
		try {
			File file = new File(filename);
			FileWriter fw = new FileWriter(file);
			fw.write("dia;S;I;R;N");

			for (int i = 0; i < matriz.length; i++) {
				fw.write("\n" + (int) matriz[i][0] + ";" + matriz[i][1] + ";" + matriz[i][2] + ";" + matriz[i][3] + ";" + matriz[i][4]);
			}
			fw.close();
			System.out.println("Ficheiro " + filename + " criado com sucesso!");

		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public static List<String[]> readFile(String filename) throws IOException {

		FileReader fr = null;
		BufferedReader br = null;
		List<String[]> result = new ArrayList<>();
		try {

			fr = new FileReader(filename);
			br = new BufferedReader((fr));
			br.readLine();
			String line;

			while ((line = br.readLine()) != null) {
				String[] values = line.split(";");
				for (int j = 0; j < values.length; j++) {
					values[j] = values[j].replace(",", ".");
				}
				result.add(values);
			}

		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			if (br != null) {
				try {
					br.close();
				} catch (IOException e) {
				}
			}
			if (fr != null) {
				try {
					fr.close();
				} catch (IOException e) {
				}
			}
		}
		return result;
	}

	public static List<Parametros> transformValues(List<String[]> values) {
		List<Parametros> parametrosList = new ArrayList<>();

		for (int i = 0; i < values.size(); i++){
			Parametros parametros = new Parametros();

			parametros.namePerson = values.get(i)[0];
			parametros.beta = Double.valueOf(values.get(i)[1]);
			parametros.gama = Double.valueOf(values.get(i)[2]);
			parametros.ro = Double.valueOf(values.get(i)[3]);
			parametros.alfa = Double.valueOf(values.get(i)[4]);

			parametrosList.add(parametros);
		}

		return parametrosList;

	}


	static class Parametros {

		public String ficheiroParametros;
		public String ficheiroResultados;
		public String metodo;
		public double stepH;
		public int populationN;
		public int daysn;
		public String ficheiroResultadosNome;
		public String namePerson;
		public double beta;
		public double gama;
		public double ro;
		public double alfa;

	}



}