#include "Mixed.cpp"

int main()
{
    try
    {
        using namespace dealii;
        using namespace MixedFEM;

        deallog.depth_console(0);
        MultithreadInfo::set_thread_limit();

        double h, pot_l2_error = 0,
                  vec_field_l2_error = 0;

        unsigned int degree		= 1;
        unsigned int l_refine = 0;

        std::ofstream prt("error_file.dat");
        for(unsigned int g_refine=2; g_refine < 7; g_refine++)
        {
            MixedPoissonProblem<2> Poisson(degree,
                                           g_refine,
                                           l_refine);

            Poisson.run_with_errors(h, pot_l2_error,
                                    vec_field_l2_error);
            prt << h << "\t"
                << pot_l2_error << "\t"
                << vec_field_l2_error << "\n";
        }
        prt.close();
    }
    catch (std::exception &exc)
    {
        std::cerr << std::endl << std::endl
                  << "----------------------------------------------------"
                  << std::endl;
        std::cerr << "Exception on processing: " << std::endl
                  << exc.what() << std::endl
                  << "Aborting!" << std::endl
                  << "----------------------------------------------------"
                  << std::endl;

        return 1;
    }
    catch (...)
    {
        std::cerr << std::endl << std::endl
                  << "----------------------------------------------------"
                  << std::endl;
        std::cerr << "Unknown exception!" << std::endl
                  << "Aborting!" << std::endl
                  << "----------------------------------------------------"
                  << std::endl;
        return 1;
    }

    return 0;
}
