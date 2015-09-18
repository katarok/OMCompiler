
#include <Core/ModelicaDefine.h>
#include <Core/Modelica.h>
#if defined(__vxworks)


#elif defined(OMC_BUILD) && !defined(RUNTIME_STATIC_LINKING)

#include <Solver/UmfPack/UmfPack.h>
#include <Solver/UmfPack/UmfPackSettings.h>

    /* OMC factory */
    using boost::extensions::factory;

BOOST_EXTENSION_TYPE_MAP_FUNCTION {
  types.get<std::map<std::string, factory<IAlgLoopSolver,IAlgLoop*, ILinSolverSettings*> > >()
    ["umfpack"].set<UmfPack>();
  types.get<std::map<std::string, factory<ILinSolverSettings> > >()
    ["umfpackSettings"].set<UmfPackSettings>();
 }
#elif defined(OMC_BUILD) && defined(RUNTIME_STATIC_LINKING)
#include <Solver/UmfPack/UmfPack.h>
#include <Solver/UmfPack/UmfPackSettings.h>

boost::shared_ptr<ILinSolverSettings> createUmfpackSettings()
{
     boost::shared_ptr<ILinSolverSettings> settings = boost::shared_ptr<ILinSolverSettings>(new UmfPackSettings());
     return settings;
}

boost::shared_ptr<IAlgLoopSolver> createUmfpackSolver(IAlgLoop* algLoop, boost::shared_ptr<ILinSolverSettings> solver_settings)
{
   boost::shared_ptr<IAlgLoopSolver> solver = boost::shared_ptr<IAlgLoopSolver>(new UmfPack(algLoop,solver_settings.get()));
   return solver;
}

#else
error "operating system not supported"
#endif



