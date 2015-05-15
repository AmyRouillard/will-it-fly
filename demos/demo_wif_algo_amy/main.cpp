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

		//std::string path = "/home/uauser/Documents/eindwerk-programmeren/will-it-fly/demos/demo_wif_algo_amy/Data/data" + std::to_string(num_lines) + ".txt";
		//std::ofstream myfile;
		//myfile.open(path);
        //myfile << "x" << "   " << "y" << "   " << "cp calc" << "   " << "cp theor pol" << "   " << "cp theor cart" << std::endl;

		std::vector<wif_core::line_2d_c> mylines = myAirfoil.get_lines();
		std::vector<wif_core::vector_2d_c> centers(num_lines);
		std::vector<double> angles(num_lines);
        std::vector<double> fouten(num_lines);
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

			//myfile << centers[i].x << "   " << centers[i].y << "   " << calculate_flow2.c_p[i] << "   " << 1 - 4 * pow(sin(angles[i]), 2) << "   " << 1 - 4 * pow((centers[i].y / radius), 2)   << std::endl;
            cp_theory[i]=1- 4 * pow((centers[i].y / radius), 2);
            fouten[i]=pow(1 - 4 * pow((centers[i].y / radius), 2)-calculate_flow2.c_p[i],2);
            x_plot[i]=i+1;

		}

        fout[q]= (std::accumulate(fouten.begin(), fouten.end(), 0))/num_lines;
		times[q]=(std::clock() - start) / (double)CLOCKS_PER_SEC ;
		som_sigma_l[q]=calculate_flow2.closed_body_check;
		num_lines_all[q]=num_lines;

		//myfile.close();

		std::vector<std::vector<double>> y_plot;
		y_plot.push_back(calculate_flow2.c_p);
		y_plot.push_back(cp_theory);

        std::vector<std::string> legende;
        legende.push_back("Numeriek bepaalde waarden");
        legende.push_back("Theoretische waarden");
        std::string title ="Cp: "+ std::to_string(num_lines)+ "panelen gebruikt";
        std::string filename= "/home/uauser/Documents/eindwerk-programmeren/will-it-fly/demos/demo_wif_algo_amy/Data/plot_cp.pdf";

        //std::shared_ptr<wif_viz::visualization_c> myRoot = wif_viz::create_visualization_root(myFlow, midpoint, midpoint);
        //myRoot->plotVectors(x_plot,y_plot, legende, filename, "Index paneel", "Cp",title);

	}

//num_lines_all,fout
//num_lines_all,times
//num_lines_all,som_sigma_l

/*    std::ofstream myfile;
    std::string path = "/home/uauser/Documents/eindwerk-programmeren/will-it-fly/demos/demo_wif_algo_amy/Data/times.txt";
    myfile.open(path);
    for (int i=0;i<times.size();i++){
    myfile << times[i]<< std::endl;
    }
    myfile.close();
*/

//////////////////////////////////////////////////////////////


	return 0;
}
