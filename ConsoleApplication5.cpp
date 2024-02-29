
#include <iostream>
#include <vector>
#include <complex>
#include <random>
#include <fstream>

using namespace std;

ofstream out;


vector<int> Arr;

class Modulator {
public:

	int bitsPerSymbol;

	vector<complex<double>> constellationPoints(int bitsPerSymbol) // Работает для всех типов
	{

		int sqrtM = pow(2, bitsPerSymbol / 2);

		vector<int> xPoints;
		for (int i = -(sqrtM - 1); i <= sqrtM - 1; i += 2) 
		{
			xPoints.push_back(i);
		}

		vector<int> yPoints;
		for (int i = sqrtM - 1; i >= -(sqrtM - 1); i -= 2) 
		{
			yPoints.push_back(i);
		}


		vector<vector<int>> x(sqrtM, xPoints);
		vector<vector<int>> y(sqrtM, yPoints);
		

		vector<int> xConstl(sqrtM * sqrtM);
		vector<int> yConstl(sqrtM * sqrtM);
		int count2 = 0;
		for (int i = 0; i < sqrtM; i++)
		{
			for (int j = 0; j < sqrtM; j++) 
			{
				xConstl[count2] = x[j][i];
				yConstl[count2] = y[i][j];
				count2++;
			}
		}


		vector<complex<double>> constellationPoints(sqrtM * sqrtM);

		for (int i = 0; i < sqrtM * sqrtM; i++)
		{
				constellationPoints[i] = complex<double>(xConstl[i], yConstl[i]);
		}
		return constellationPoints;
	}

	vector<int> gray(int bitsPerSymbol, vector<complex<double>> Points){
		vector<int> binMapping;
		for (int i = 0; i < pow(2, bitsPerSymbol); i++)
		{
			binMapping.push_back(i);
		}

		int nBitsBy2 = bitsPerSymbol >> 1;
		vector<int> symbolI;
		vector<int> symbolQ;

		for (int i = 0; i < binMapping.size(); i++)
		{
			symbolI.push_back(binMapping[i] >> nBitsBy2); 
		}

		for (int i = 0; i < binMapping.size(); i++)
		{
			int a = pow(2, bitsPerSymbol) - 1;
			int b = a >> nBitsBy2;
			symbolQ.push_back(binMapping[i] & b); 
		}

		int i = 1;
		vector<int> tmpI;
		vector<int> tmpQ;
		vector<int> symboltmpI = symbolI;
		vector<int> symboltmpQ = symbolQ;

		while (i < nBitsBy2)
		{
			for (int j = 0; j < symbolI.size(); j++)
			{
				tmpI.push_back(symboltmpI[j] >> i);
				tmpQ.push_back(symboltmpQ[j] >> i);

				symboltmpI[j] = symboltmpI[j] ^ tmpI[j];
				symboltmpQ[j] = symboltmpQ[j] ^ tmpQ[j];

			}
			i = i + 1;
			tmpI.clear();
			tmpQ.clear();
		}

		vector<int> SymbolIndex;
		for (int i = 0; i < binMapping.size(); i++)
		{
			int a = symboltmpI[i] << nBitsBy2;
			int b = a + symboltmpQ[i];
			SymbolIndex.push_back(b);
		}
		return SymbolIndex;
	}


	vector<complex<double>> modulate(const vector<int>& data) 
	{
		vector<int> powOf2;

		for (int i = bitsPerSymbol - 1; i >= 0; i--)
		{
			int temp = pow(2, i);
			powOf2.push_back(temp);
		}


		vector<vector<int>> bits(data.size() / bitsPerSymbol, vector <int>(bitsPerSymbol)); 
		int count = 0;
		for (int i = 0; i < data.size() / bitsPerSymbol; i++)
		{
			for (int j = 0; j < bitsPerSymbol; j++)
			{
				bits[i][j] = data[count];
				count++;
			}
		}

		vector<vector<int>> bitInputMat(data.size() / bitsPerSymbol, vector <int>(bitsPerSymbol)); 
		for (int j = 0; j < bitsPerSymbol; j++)
			for (int i = 0; i < data.size() / bitsPerSymbol; i++)
			{
				{
					bitInputMat[i][j] = bits[i][j] * powOf2[j];
				}
			}

		vector<int> intOutput(data.size() / bitsPerSymbol);

		int count3 = 0;
		int sum;
		for (int i = 0; i < data.size() / bitsPerSymbol; i++) 
		{
			sum = 0;
			for (int j = 0; j < bitsPerSymbol; j++)
			{ 
				sum = sum + bitInputMat[i][j];
				intOutput[count3] = sum;
			}
			count3++;
		}
	
		vector<complex<double>> Points = constellationPoints(bitsPerSymbol);
		vector<int> SymbolIndex = gray(bitsPerSymbol, Points);

		vector<double> modulatedSignal1(intOutput.size());
		for (int i = 0; i < modulatedSignal1.size(); i++) 
		{
			modulatedSignal1[i] = SymbolIndex[intOutput[i]];
		}

			vector<complex<double>> modulatedSignal(intOutput.size());

			for (int i = 0; i < modulatedSignal.size(); i++)
			{
				modulatedSignal[i] = Points[modulatedSignal1[i]];
			}

		return modulatedSignal;
	}
};


class Channel: public Modulator {
public:
	double EbNo;

public:
	vector<complex<double>> AddNoise(vector<complex<double>> dataInput, double EbNo)
	{
		vector<complex<double>> Points = Modulator::constellationPoints(Modulator::bitsPerSymbol);
		
		double temp = 0;

		for (int i = 0; i < Points.size(); i++)
		{
			temp = pow((abs(Points[i])), 2) + temp;
		}
		double PS = (double) temp / Points.size();
		double Pbd = (double) PS / Modulator::bitsPerSymbol;

		double Sigma = sqrt(Pbd* pow(10, -EbNo / 10) / 2);

		normal_distribution<> distr(0, 1);
		random_device gen;
		vector<double> ndRe(dataInput.size());
		vector<double> ndIm(dataInput.size());

		for (int i = 0; i < dataInput.size(); i++) 
		{
			ndRe[i] = distr(gen);
			ndIm[i] = distr(gen);
		}
		vector<complex<double>> NoiseArr(dataInput.size());
		vector<complex<double>> Noise(dataInput.size());
		for (int i = 0; i < dataInput.size(); i++)
		{
			NoiseArr[i] = complex<double>(ndRe[i], ndIm[i]);
			Noise[i] = NoiseArr[i] * Sigma;
		}
		
		vector<complex<double>> SigAndNoise;
		for (int i = 0; i < dataInput.size(); i++)
		{
			SigAndNoise.push_back(dataInput[i] + Noise[i]);
		}
		//vector<complex<double>> SigAndNoise = dataInput;
		return SigAndNoise;
	}
};

class Demodulator: public Modulator{
	
public:
	
	vector<int> demodulate(vector<complex<double>> dataInput)
	{
		int sqrtM = pow(2, bitsPerSymbol / 2);

		vector<int> rIdx(dataInput.size());
		for (int i = 0; i < dataInput.size(); i++) {
			rIdx[i] = round(((dataInput[i].real() + (sqrtM - 1)) / 2));
			if (rIdx[i] < 0) {
				rIdx[i] = 0;
			}
			if (rIdx[i] > (sqrtM - 1)) {
				rIdx[i] = sqrtM - 1;
			}
		}

		
		vector<int> iIdx(dataInput.size());
		for (int i = 0; i < dataInput.size(); i++) {
			iIdx[i] = round(((dataInput[i].imag() + (sqrtM - 1)) / 2));
			if (iIdx[i] < 0) {
				iIdx[i] = 0;
			}
			if (iIdx[i] > (sqrtM - 1)) {
				iIdx[i] = sqrtM - 1;
			}
		}

		vector<int> z;
		for (int i = 0; i < dataInput.size(); i++) 
		{
			int temp1 = sqrtM - iIdx[i] - 1;
			int temp2 = sqrtM * rIdx[i];
			z.push_back(temp1 + temp2);
		}

		vector<complex<double>> Points = Modulator::constellationPoints(Modulator::bitsPerSymbol);
		
		vector<int> SymbolIndex = Modulator::gray(Modulator::bitsPerSymbol, Points);

		vector<int> mapping(SymbolIndex.size());
		for (int i = 0; i < pow(2, Modulator::bitsPerSymbol); i++)
		{
			mapping[SymbolIndex[i]] = i;
		}


		vector<double> modulatedSignal1(z.size());
		for (int i = 0; i < modulatedSignal1.size(); i++)
		{
			modulatedSignal1[i] = mapping[z[i]];
		}


		vector<int> bitOutput(modulatedSignal1.size()*bitsPerSymbol);

		vector<double> powOf2;

		for (int i = -bitsPerSymbol + 1; i <= 0; i++)
		{
			double temp = pow(2, i);
			powOf2.push_back(temp);
		}

		vector<vector<double>> bitOutputTmp(modulatedSignal1.size(), vector <double>(bitsPerSymbol));
		
		for (int i = 0; i < bitsPerSymbol; i++)
		{
			int count = 0;
			for (int j = 0; j < modulatedSignal1.size(); j++)
			{
				bitOutputTmp[j][i] = modulatedSignal1[count];
				count++;
			}
		}

		vector<vector<double>> bitOutputTmp1(modulatedSignal1.size(), vector <double>(bitsPerSymbol));

		for (int i = 0; i < bitsPerSymbol; i++)
		{
			for (int j = 0; j < modulatedSignal1.size(); j++)
			{
				bitOutputTmp1[j][i] = floor(bitOutputTmp[j][i]* powOf2[i]);
			}
		}

		vector<vector<int>> bitOutputTmpInt(bitOutputTmp1.size(), vector <int>(bitsPerSymbol));

		for (int i = 0; i < bitOutputTmp1.size(); i++)
		{
			for (int j = 0; j < bitsPerSymbol; j++)
			{
				int temp = bitOutputTmp1[i][j];
				bitOutputTmpInt[i][j] = temp % 2;
			}
		}

		int count2 = 0;
		for (int i = 0; i < bitOutputTmpInt.size(); i++)
		{
			for (int j = 0; j < bitsPerSymbol; j++)
			{
				bitOutput[count2] = bitOutputTmpInt[i][j];
				count2++;
			}
		}

		return bitOutput;
	}
};

int main()
{
	out.open("BER_results.txt");


	int lenbitsstream = 38400;
	for (int i = 0; i < lenbitsstream; i++)
	{
		int temp = rand() % (2);
		Arr.push_back(temp);
	}

	vector<double> EbNo;
	
	int minSNR = 0;
	int stepSNR = 1;
	int maxSNR = 10;

	int count = 0;
	for (int j = minSNR; j <= maxSNR; j = j + stepSNR) 
	{
		EbNo.push_back(j);
		count++;
	}

	vector<int> ModOrder = {2, 4, 6, 8}; //QAM4, QAM16, QAM64, QAM256
	Modulator Modqam4;
	Channel channel;
	Demodulator Demodqam4;

	
	vector<vector<double>> BER(EbNo.size(), vector <double>(ModOrder.size()));
	vector<double> BER_count(EbNo.size());

	//int c = 0;
	for (int jj = 0; jj < ModOrder.size(); jj++)
	{
		Modqam4.bitsPerSymbol = ModOrder[jj];
		channel.bitsPerSymbol = ModOrder[jj];
		Demodqam4.bitsPerSymbol = ModOrder[jj];


		for (int j = 0; j < EbNo.size(); j++)
		{

			channel.EbNo = EbNo[j];

			vector<complex<double>> ModSig = Modqam4.modulate(Arr);
			vector<complex<double>> SigAndNoise = channel.AddNoise(ModSig, channel.EbNo);
			vector<int> DeModSig = Demodqam4.demodulate(SigAndNoise);


			int count = 0;

			for (int i = 0; i < Arr.size(); i++)
			{
				if (Arr[i] == DeModSig[i])
				{
					count = count;
				}
				else
				{
					count++;
				}
			}
			int Len = Arr.size();
			double err = (double)count / Len;
			BER_count[j] = err;
		}
		for (int ii = 0; ii < EbNo.size(); ii++) 
		{
			BER[ii][jj] = BER_count[ii];
			out << BER[ii][jj];
			out << endl;
		}
	}
	out.close();
}