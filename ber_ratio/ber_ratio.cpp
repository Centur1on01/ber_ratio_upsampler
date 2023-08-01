
#include <iostream>
#include <complex>
#include <vector>
#include <ctime>
#include <cmath>
#include <random>
#include <string>
#include <fstream>

using namespace std;

class gen {
    double real, imag;
    complex<double> result;

public:

    gen() {}

    gen(double& Ps) {
        real = sqrt(Ps / 2.);
        imag = sqrt(Ps / 2.);
        result.real(real);
        result.imag(imag);
        srand(time(0));
    }

    complex<double> next() {
        return (rand() % 2 == 0 ? result : -result);
    }

};

class awgn {
    default_random_engine generator;
    normal_distribution<double> distr;
    double Pn;
public:

    awgn(const double& m, const double& s)
    {
        Pn = s;
        distr = normal_distribution<double>(m, 1);
    }

    double next() {
        return std::sqrt(Pn / 2) * distr(generator);
    }

};

class upsampling {
    double n;

public:

    upsampling() {}

    upsampling(double& n) {
        this->n = n;
    }

    complex<double> next(complex<double>& signal, awgn& noise) {
        signal.real(signal.real() + noise.next() / n);
        signal.imag(signal.imag() + noise.next() / n);
        return signal;
    }

};



class decisionmaker {
    complex<double> result;
    double Ps;

public:

    decisionmaker() {}

    decisionmaker(double Ps) {
        this->Ps = Ps;
    }

    complex<double> decide(complex<double> signal) {
        if (signal.imag() >= -signal.real()) {
            result.real(sqrt(Ps / 2.));
            result.imag(result.real());
            return result;
        }
        else {
            result.real(-sqrt(Ps / 2.));
            result.imag(result.real());
            return result;
        }
    }

};


int main()
{
    double Ps = 2.;
    gen generator(Ps);
    ofstream signalout("out1.txt");
    ofstream signalnoisedout("out1n.txt");
    ofstream signaldecisionmade("out1d.txt");
    ofstream signalberratio("out2.txt");
    ofstream signalupscaledout("outup1.txt");
    complex<double> signal;
    complex<double> signaln;
    complex<double> signald;
    complex<double> signalupn;
    complex<double> signalupd;
    double signalupn_real = 0;
    double signalupn_imag = 0;
    double ampsum = 0;
    //double ampsum_upn = 0;
    double correct = 0;
    double incorrect = 0;
    double PsdB = 10. * log10(Ps);
    double SNR = 0;
    double PndB = PsdB - SNR;
    double BER_t;
    double N = 100000;
    //awgn noise(0, PsdB - SNR);
    double upsampling_n = 1; // 1/2/4/8
    vector<complex<double>> buff;
    upsampling upsampler(upsampling_n);
    decisionmaker decision(Ps);
    int it = 0;
    for (int snr = 0; snr <= 10; ++snr) {
        PndB = PsdB - snr;
        double Pn = pow(10., PndB / 10.);
        awgn noise(0, Pn);
        ampsum = 0;
        //ampsum_upn = 0;
        correct = 0;
        incorrect = 0;
        for (int i = 0; i < N; ++i) {
            signal = generator.next();
            signaln = signal;
            for (int j = 0; j < upsampling_n; ++j) {
                signalupn = upsampler.next(signaln, noise);
                buff.push_back(signalupn);
                signalupn_real += buff[j].real(); // сумма всех повторений сигнала
                signalupn_imag += buff[j].imag();
                if (j == upsampling_n - 1) {
                    signalupd.real(signalupn_real / upsampling_n); // среднее арифметическое
                    signalupd.imag(signalupn_imag / upsampling_n);
                    signalupd = decision.decide(signalupd);
                    if (signalupd == signal) {
                        ++correct;
                    }
                    //signalupscaledout << i << "\t" << signalupd.real() << (signalupd.imag() < 0 ? "" : "+") << signalupd.imag() << "i" << endl;
                    for (int k = 0; k < upsampling_n; ++k) {
                        ampsum += pow(signalupd.real() - buff[k].real(), 2) + pow(signalupd.imag() - buff[k].imag(), 2); // вектор разности считается для каждого k повторения сигнала
                    }
                    //ampsum_upn += pow(signalupd.real() - signalupn_real, 2) + pow(signalupd.imag() - signalupn_imag, 2);
                    buff.clear();
                    signalupd.real(0);
                    signalupd.imag(0);
                    signalupn_real = 0;
                    signalupn_imag = 0;
                }

            }

            /*signaln.real(signaln.real() + noise.next());
            signaln.imag(signaln.imag() + noise.next());
            signald = decision.decide(signaln);
            if (it == 0) {
                signalout << i << "\t" << signal.real() << (signal.imag() < 0 ? "" : "+") << signal.imag() << "i" << endl;
                signalnoisedout << i << "\t" << signaln.real() << (signaln.imag() < 0 ? "" : "+") << signaln.imag() << "i" << endl;
                signaldecisionmade << i << "\t" << signald.real() << (signald.imag() < 0 ? "" : "+") << signald.imag() << "i" << endl;
            }
            ampsum += pow(signald.real() - signaln.real(), 2) + pow(signald.imag() - signaln.imag(), 2);
            if (signald == signal)
                ++correct;*/
        }
        BER_t = 0.5 * erfc(sqrt(pow(10., snr / 10.)));
        signalberratio << snr << "\t" << (N - correct) / N << "\t" << BER_t << endl;
        ++it;
        cout << "Ps = " << Ps << " | Ps dB = " << PsdB << " | SNR dB = " << snr << endl;
        cout << "Pn = " << pow(10., PndB / 10.) << " | Pn dB = " << PndB << endl;
        cout << "Pn calculated = " << ampsum / N << " | dB = " << 10. * log10(ampsum / N) << endl;
        cout << "correct = " << correct << " | incorrect = " << N - correct << " | BER = " << (N - correct) / N << endl << endl;
    }

    signalout.close();
    signalnoisedout.close();
    signaldecisionmade.close();
    system("pause");
}
