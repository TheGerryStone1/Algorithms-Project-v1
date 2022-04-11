#include <iostream>
#include <string>
#include <vector>
#include <fstream>

using namespace std;

ofstream outfile;

vector<int> lpsvFunction(string p)
{
    int n = p.length(), i = 1, j = 0;
    vector<int> lpsv(n, 0);
    
    while(i < n)
    {
        if(p[i] == p[j])
        {
            lpsv[i] = j + 1;
            i++;
            j++;
        }

        else
        {
            if(j == 0)
            {
                lpsv[i] = 0;
                i++;
            }
            else j = lpsv[j - 1];
        }
    }

    return lpsv;
}

vector<int> kmpAlgorithm(string text, string pattern)
{
    vector<int> posMatch;
    vector<int> lpsv = lpsvFunction(pattern);
    int i = 0, j = 0, n = text.length(), m = pattern.length();

    while(i < n)
    {
        if(text[i] == pattern[j])
        {
            i++;
            j++;
            if(j == m)
            {
                posMatch.push_back(i - m);
                j = lpsv[j - 1];
            }
        }

        else
        {
            if(j == 0) i++;
            else j = lpsv[j - 1];
        }
    }

    return posMatch;
}

int lcsArr[1001][1001];

string longestCommonSubstring(string transmission1, string transmission2, int size1, int size2)
{
    int longest = 0, finalIndex;
    string word  = "";
    
    lcsArr[0][0] = 0;
    for(int column = 1; column < size2 + 1; column++) lcsArr[0][column] = 0;
    for(int row = 1; row < size1 + 1; row++) lcsArr[row][0] = 0;

    for (int i = 1; i < size1 + 1; i++)
    {
        for (int j = 1; j < size2 + 1; j++)
        {
            if(transmission1[i - 1] == transmission2[j - 1]) 
            {
                lcsArr[i][j] = lcsArr[i - 1][j - 1] + 1;
                longest = max(longest, lcsArr[i][j]);
            }

            else lcsArr[i][j] = 0;
        }
    }

    for (int i = 1; i < size1 + 1; i++)
    {
        for (int j = 1; j < size2 + 1; j++)
        {
            if(lcsArr[i][j] == longest) 
            {
                finalIndex = i - 1;
            }
        }
    }

    for(int k = finalIndex; k > finalIndex - longest; k--)
    {
        word = transmission1[k] + word;
    }

    return word;
}

string manacher(string texto)
{
		string T = "";
		int n = texto.length();

		for (int i=0; i < n; i++)
		{
			T += "|";
			T += texto[i];
		}

		T += "|";
		n = T.length();
		vector<int> L(n);
		L[0] = 0;
		L[1] = 1;
		int maxLong = 0, maxCentro = 0;
		bool exp;

		for (int c = 1, li = 0, ri = 2; ri < n; ri++)
		{
			li = c - (ri - c);
			exp = false;

			if (c - maxLong <= ri && ri >= c + maxLong){
				if (L[li] < (c + L[c] - ri)){
					L[ri] = L[li];
				}

				else if(L[li] == (c + L[c]) - ri && (c + L[c]) == 2 * n)
				{
					L[ri] = L[li];
				}

				else if(L[li] == (c + L[c]) - ri && (c + L[c]) < 2 * n)
				{
					L[ri] = L[li];
					exp = true;
				}

				else if(L[li] > (c + L[c]) - ri && (c + L[c]) < 2 * n)
				{
					L[ri] = (c+L[c]) - ri;
					exp = true;
				}
			}

			else
			{
				L[ri] = 0;
				exp = true;
			}

			if (exp)
			{
				while((ri + L[ri] < n) && (ri - L[ri] > 0) && (T[ri - L[ri] - 1] == T[ri + L[ri] + 1]))
				{
					L[ri]++;
				}
			}

			c = ri;
			if (L[ri] > maxLong)
			{
				maxLong = L[ri];
				maxCentro = ri;
			}
		}

		int inicio = (maxCentro - maxLong) / 2;
		string salida = "";

		for (int i = inicio; i < (inicio + maxLong); i++)
		{
			salida += texto[i];
		}

        outfile << "Posicion: " << inicio << endl;
		return salida;
}

int main()
{
    string x;
    vector<string> transmissionsContent;
    vector<string> transmissionsNames;
    vector<string> codes;
    fstream newFile;

    outfile.open("checking.txt");

    //Read codes
    newFile.open("mcode.txt");
    while(newFile >> x)
    {
        codes.push_back(x);    
    }
    newFile.close();

    //Read transmissions
    string tempTransmissionName, line, tempTransmissionContent;
    for(int i = 1; i <= 3; i++)
    {
        tempTransmissionName = "transmission" + to_string(i) + ".txt";
        newFile.open(tempTransmissionName);
        while(getline(newFile, line))
        {    
            tempTransmissionContent += line;
        }
        newFile.close();
        transmissionsNames.push_back(tempTransmissionName);
        transmissionsContent.push_back(tempTransmissionContent);

        tempTransmissionContent = "";
    }

    //FIND ALL ITERATIONS OF MALICIOUS CODE || FIND ALL ITERATIONS OF MALICIOUS CODE || FIND ALL ITERATIONS OF MALICIOUS CODE
    //FIND ALL ITERATIONS OF MALICIOUS CODE || FIND ALL ITERATIONS OF MALICIOUS CODE || FIND ALL ITERATIONS OF MALICIOUS CODE
    //FIND ALL ITERATIONS OF MALICIOUS CODE || FIND ALL ITERATIONS OF MALICIOUS CODE || FIND ALL ITERATIONS OF MALICIOUS CODE

    vector<int> posMatch;

    for(int j = 0; j < codes.size(); j++)
    {
        outfile << "Codigo: " << codes[j] << endl;

        for(int i = 0; i < transmissionsContent.size(); i++)
        {
            posMatch = kmpAlgorithm(transmissionsContent[i], codes[j]);

            outfile << transmissionsNames[i] << " ====> ";

            if(posMatch.size() == 0) outfile << "No match";
            else
            {
                outfile << posMatch.size();

                if(posMatch.size() == 1) outfile << " vez." << endl;
                else outfile << " veces." << endl;

                for(int i = 0; i < posMatch.size(); i++)
                {
                    outfile << posMatch[i];
                    if(i != posMatch.size() - 1) outfile << ", ";
                }
            }
            outfile << endl;
        }

        if(j != codes.size() - 1) outfile << "-------------------------------------------------------------------------" << endl;
        else outfile << endl;
    }

    outfile << "=========================================================================" << endl;
    outfile << "=========================================================================" << endl;
    outfile << endl;

    //LONGEST PALINDROME || LONGEST PALINDROME || LONGEST PALINDROME || LONGEST PALINDROME
    //LONGEST PALINDROME || LONGEST PALINDROME || LONGEST PALINDROME || LONGEST PALINDROME
    //LONGEST PALINDROME || LONGEST PALINDROME || LONGEST PALINDROME || LONGEST PALINDROME

    outfile << "Palindromo mas grande: " << endl;

    for(int i = 0; i < transmissionsContent.size(); i++)
    {
        outfile << transmissionsNames[i] << " ====> " << manacher(transmissionsContent[i]) << endl;
    }

    outfile << endl;
    outfile << "=========================================================================" << endl;
    outfile << "=========================================================================" << endl;
    outfile << endl;

    //LONGEST COMMON SUBSTRING || LONGEST COMMON SUBSTRING || LONGEST COMMON SUBSTRING || LONGEST COMMON SUBSTRING
    //LONGEST COMMON SUBSTRING || LONGEST COMMON SUBSTRING || LONGEST COMMON SUBSTRING || LONGEST COMMON SUBSTRING
    //LONGEST COMMON SUBSTRING || LONGEST COMMON SUBSTRING || LONGEST COMMON SUBSTRING || LONGEST COMMON SUBSTRING

    string transmission1 = transmissionsContent[0], transmission2 = transmissionsContent[1], transmission3 = transmissionsContent[2];
    int size1 = transmission1.length(), size2 = transmission2.length(), size3 = transmission3.length();
    string lcs1 = longestCommonSubstring(transmission1, transmission2, size1, size2);
    string lcs2 = longestCommonSubstring(transmission1, transmission3, size1, size3);
    string lcs3 = longestCommonSubstring(transmission2, transmission3, size2, size3);
    
    if(lcs1.size() >= lcs2.size() && lcs1.size() >= lcs3.size()) outfile << "Substring mas largo: " << lcs1 << endl;
    else if (lcs2.size() >= lcs1.size() && lcs2.size() >= lcs3.size()) outfile << "Substring mas largo: " << lcs2 << endl;
    else outfile << "Substring mas largo: " << lcs3 << endl;

    outfile.close();
}