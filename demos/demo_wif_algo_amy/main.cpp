#include <iostream>
#include <cmath>
#include <fstream>

#include <wif_core/wif_core.hpp>
#include <wif_algo/wif_algo.hpp>
#include <wif_viz/wif_viz.hpp>

#include <ctime>
int main()
{

	double pi = 3.1415;
	double radius = 2;


	unsigned int n=20;
	std::clock_t start;
	std::vector<double> times(n);
	std::vector<double> fout(n);
	std::vector<double> som_sigma_l(n);
	std::vector<double> num_lines_all(n);

	for(int q = 1; q < n+1; q++)
	{
		start = std::clock();

		unsigned int num_lines = 10 * q;
		wif_core::vector_2d_c midpoint(0, 0);
		wif_core::airfoil_c myAirfoil(midpoint, radius, num_lines);
		std::shared_ptr<wif_core::uniform_flow_c> myFlow = std::make_shared<wif_core::uniform_flow_c>(0, 1);
		bool Kutta = 0;
		wif_algo::calculation_results_c calculate_flow2 = wif_algo::calculate_flow(myAirfoil, myFlow, Kutta);
		double gamma = 0;

		std::vector<wif_core::line_2d_c> mylines = myAirfoil.get_lines();
		std::vector<wif_core::vector_2d_c> centers(num_lines);
		std::vector<double> angles(num_lines);
	    double fout_accum=0;
	    std::vector<double> x_plot(num_lines);
	    std::vector<double> cp_theory(num_lines);

		for(unsigned int i = 0; i < num_lines; i++)
		{
			wif_core::line_2d_c temp_line = mylines[i];
			centers[i] = temp_line.get_center_point();

			if(centers[i].y > 0)
			{
				angles[i] = atan2(centers[i].y, centers[i].x);
			}
			else
			{
				angles[i] = atan2(centers[i].y, centers[i].x) + 2 * pi;
			}

			cp_theory[i]=1- 4 * pow((centers[i].y / radius), 2);
	        fout_accum+=pow(1 - 4 * pow((centers[i].y / radius), 2)-calculate_flow2.c_p[i],2);
	        x_plot[i]=i+1;

		}

	    fout[q-1]= log10((fout_accum/num_lines));
		times[q-1]=(std::clock() - start) / (double)CLOCKS_PER_SEC ;
		som_sigma_l[q-1]=(calculate_flow2.closed_body_check)*1000;
		num_lines_all[q-1]=num_lines;

	    std::cout<< fout[q-1] <<std::endl;

		std::vector<std::vector<double>> y_plot;
		y_plot.push_back(calculate_flow2.c_p);
		y_plot.push_back(cp_theory);

	    std::vector<std::string> legende_cp;
	    legende_cp.push_back("Numeriek bepaalde waarden");
	    legende_cp.push_back("Theoretische waarden");
	    std::string title_cp ="Cp: "+ std::to_string(num_lines)+ "panelen gebruikt";
	    std::string filename_cp= "/home/uauser/Documents/eindwerk-programmeren/will-it-fly/demos/demo_wif_algo_amy/Data/plot_cp "+std::to_string(num_lines)+".png";

	    std::shared_ptr<wif_viz::visualization_c> myRoot_1 = wif_viz::create_visualization_root(myFlow, midpoint, midpoint);
	    myRoot_1->plotVectors(y_plot,x_plot, legende_cp, filename_cp, "Index paneel", "Cp",title_cp);

	}

    //Dummy myFlow en midpoint:
    wif_core::vector_2d_c midpoint(0, 0);
    wif_core::airfoil_c myAirfoil(midpoint, 2, 10);
    std::shared_ptr<wif_core::uniform_flow_c> myFlow = std::make_shared<wif_core::uniform_flow_c>(0, 1);

    //Plots:
    std::vector<std::string> legende_error;
    legende_error.push_back("Fout");
    std::string title_error ="Fout op cp als functie van het aantal panelen";
    std::string filename_error= "/home/uauser/Documents/eindwerk-programmeren/will-it-fly/demos/demo_wif_algo_amy/Data/fout op cp.png";
    std::shared_ptr<wif_viz::visualization_c> myRoot_2 = wif_viz::create_visualization_root(myFlow, midpoint, midpoint);
    std::vector<std::vector<double>> y_plot_fout;
    y_plot_fout.push_back(fout);
    myRoot_2 ->plotVectors(y_plot_fout,num_lines_all, legende_error, filename_error, "Aantal panelen", "log10(Fout op Cp)",title_error);

    std::vector<std::string> legende_time;
    legende_time.push_back("Tijdsduur");
    std::string title_time ="Tijdsduur als functie van het aantal panelen";
    std::string filename_time= "/home/uauser/Documents/eindwerk-programmeren/will-it-fly/demos/demo_wif_algo_amy/Data/tijd cp.png";
    std::shared_ptr<wif_viz::visualization_c> myRoot_3 = wif_viz::create_visualization_root(myFlow, midpoint, midpoint);
    std::vector<std::vector<double>> y_plot_time;
    y_plot_time.push_back(times);
    myRoot_3 ->plotVectors(y_plot_time,num_lines_all,legende_time, filename_time, "Aantal panelen", "Tijd (s)",title_time);

    std::vector<std::string> legende_sigma;
    legende_sigma.push_back("Som sigma * l");
    std::string title_sigma ="Som sigma * l als functie van het aantal panelen";
    std::string filename_sigma= "/home/uauser/Documents/eindwerk-programmeren/will-it-fly/demos/demo_wif_algo_amy/Data/som sigma l cp.png";
    std::shared_ptr<wif_viz::visualization_c> myRoot_4 = wif_viz::create_visualization_root(myFlow, midpoint, midpoint);
    std::vector<std::vector<double>> y_plot_sigma;
    y_plot_sigma.push_back(som_sigma_l);
    myRoot_4 ->plotVectors(y_plot_sigma,num_lines_all,legende_sigma, filename_sigma, "Aantal panelen", "Som sigma *l (10^3)",title_sigma);



    ///////////////////////////////////////////////////

    double radius_airfoils = 1;
    unsigned int num_lines_airfoils = 200;
    bool Kutta_airfoils = 1;
    std::string path_airfoils="/home/uauser/Documents/eindwerk-programmeren/will-it-fly/demos/demo_wif_algo_amy/";
    wif_core::vector_2d_c midpoint_airfoils(0, 0);
    //std::cout << "radius = " << radius_airfoils << std::endl;

    std::vector<double> c_l_1;
    std::vector<double> c_l_2;
    std::vector<double> c_l_3;
    std::vector<double> c_l_4;

    std::vector<double> gamma_1;
    std::vector<double> gamma_2;
    std::vector<double> gamma_3;
    std::vector<double> gamma_4;

     std::vector<double> angles;

    for (int angle_count=0; angle_count<19; angle_count++){

        double angle_in_degrees = angle_count*10-90;
        angles.push_back(angle_in_degrees);

        std::shared_ptr<wif_core::uniform_flow_c> UFlow = std::make_shared<wif_core::uniform_flow_c>((angle_in_degrees / 180) * pi, 1);
        std::cout << "angle of attack in degrees = " << angle_in_degrees << std::endl;

        //Verschillende airfoils:
        //A18 http://airfoiltools.com/airfoil/details?airfoil=a18-il
        std::cout<<"A18 (original)"<<std::endl;
        wif_core::airfoil_c Airfoil_1(path_airfoils+"A18.dat");
        Airfoil_1 = Airfoil_1.get_circle_projection(num_lines_airfoils, midpoint_airfoils, radius_airfoils, 0.0001);
        Airfoil_1 = Airfoil_1.closed_merge();
        wif_algo::calculation_results_c results_Airfoil_1 = wif_algo::calculate_flow(Airfoil_1, UFlow, Kutta_airfoils);
        std::cout<<"c_l = "<<results_Airfoil_1.c_l << " gamma = "<< results_Airfoil_1.gamma <<std::endl;
        c_l_1.push_back(results_Airfoil_1.c_l);
        gamma_1.push_back(results_Airfoil_1.gamma);

        //NACA 0012 http://airfoiltools.com/airfoil/details?airfoil=n0012-il
        std::cout<<"NACA 0012"<<std::endl;
        wif_core::airfoil_c Airfoil_2(path_airfoils+"NACA 0012.dat");
        Airfoil_2 = Airfoil_2.get_circle_projection(num_lines_airfoils, midpoint_airfoils, radius_airfoils, 0.0001);
        Airfoil_2 = Airfoil_2.closed_merge();
        wif_algo::calculation_results_c results_Airfoil_2 = wif_algo::calculate_flow(Airfoil_2, UFlow, Kutta_airfoils);
        std::cout<<"c_l = "<<results_Airfoil_2.c_l<< " gamma = "<< results_Airfoil_2.gamma<<std::endl;
        c_l_2.push_back(results_Airfoil_2.c_l);
        gamma_2.push_back(results_Airfoil_2.gamma);

        //DAVIS BASIC B-24 WING http://airfoiltools.com/plotter/index?airfoil=davis_corrected-il
        std::cout<<"DAVIS BASIC B-24 WING"<<std::endl;
        wif_core::airfoil_c Airfoil_3(path_airfoils+"DAVIS BASIC B-24 WING.dat");
        Airfoil_3 = Airfoil_3.get_circle_projection(num_lines_airfoils, midpoint_airfoils, radius_airfoils, 0.0001);
        Airfoil_3 = Airfoil_3.closed_merge();
        wif_algo::calculation_results_c results_Airfoil_3 = wif_algo::calculate_flow(Airfoil_3, UFlow, Kutta_airfoils);
        std::cout<<"c_l = "<<results_Airfoil_3.c_l<< " gamma = "<< results_Airfoil_3.gamma<<std::endl;
        c_l_3.push_back(results_Airfoil_3.c_l);
        gamma_3.push_back(results_Airfoil_3.gamma);

        //PMC19 Smoothed http://airfoiltools.com/airfoil/lednicerdatfile?airfoil=pmc19sm-il
        std::cout<<"PMC19 Smoothed"<<std::endl;
        wif_core::airfoil_c Airfoil_4(path_airfoils+"PMC19 Smoothed.dat");
        Airfoil_4 = Airfoil_4.get_circle_projection(num_lines_airfoils, midpoint_airfoils, radius_airfoils, 0.0001);
        Airfoil_4 = Airfoil_4.closed_merge();
        wif_algo::calculation_results_c results_Airfoil_4 = wif_algo::calculate_flow(Airfoil_4, UFlow, Kutta_airfoils);
        std::cout<<"c_l = "<<results_Airfoil_4.c_l << " gamma = "<< results_Airfoil_4.gamma <<std::endl;
        c_l_4.push_back(results_Airfoil_4.c_l);
        gamma_4.push_back(results_Airfoil_4.gamma);

}

    std::vector<std::vector<double>> c_l;
    c_l.push_back(c_l_1);
    c_l.push_back(c_l_2);
    c_l.push_back(c_l_3);
    c_l.push_back(c_l_4);

    std::vector<std::vector<double>> gamma;
    gamma.push_back(gamma_1);
    gamma.push_back(gamma_2);
    gamma.push_back(gamma_3);
    gamma.push_back(gamma_4);

    //Plot:
    //Dummy Uflow
    std::shared_ptr<wif_core::uniform_flow_c> UFlow = std::make_shared<wif_core::uniform_flow_c>(0, 1);
    //

    std::vector<std::string> legende_c_l;
    legende_c_l.push_back("A18");
    legende_c_l.push_back("NACA 0012");
    legende_c_l.push_back("DAVIS BASIC B-24 WING");
    legende_c_l.push_back("PMC19 Smoothed");

    std::string title_c_l ="c_l bij verschillende 'angles of attack' voor verschillende airfoil";
    std::string filename_c_l= "/home/uauser/Documents/eindwerk-programmeren/will-it-fly/demos/demo_wif_algo_amy/Data/c_l voor verschillende airfoils.png";
    std::shared_ptr<wif_viz::visualization_c> myRoot_c_l = wif_viz::create_visualization_root(UFlow, midpoint_airfoils, midpoint_airfoils);
    myRoot_c_l ->plotVectors(c_l,angles,legende_c_l, filename_c_l, "Hoek (graden)", "Cl, lift coefficent",title_c_l);

    std::string title_gamma ="gamma bij verschillende 'angles of attack' voor verschillende airfoil";
    std::string filename_gamma= "/home/uauser/Documents/eindwerk-programmeren/will-it-fly/demos/demo_wif_algo_amy/Data/gamma voor verschillende airfoils.png";
    std::shared_ptr<wif_viz::visualization_c> myRoot_gamma = wif_viz::create_visualization_root(UFlow, midpoint_airfoils, midpoint_airfoils);
    myRoot_gamma ->plotVectors(gamma,angles,legende_c_l, filename_gamma, "Hoek (graden)", "Gamma",title_gamma);

	return 0;
}
