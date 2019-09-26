﻿using System.Collections.Generic;
using MGroup.FEM.Elements;
using MGroup.FEM.Entities;
using MGroup.FEM.Interpolation;
using MGroup.FEM.Interpolation.GaussPointExtrapolation;
using MGroup.LinearAlgebra.Matrices;
using MGroup.Materials;
using MGroup.MSolve.Discretization.Integration;
using MGroup.MSolve.Discretization.Integration.Quadratures;
using MGroup.MSolve.Discretization.Mesh;
using TriangleNet.Tools;
using Xunit;

//TODO: Add tests for wrong node orders, too distorted shapes, etc.
//TODO: Add tests presented in https://www.colorado.edu/engineering/CAS/courses.d/IFEM.d/IFEM.Ch24.d/IFEM.Ch24.pdf
namespace MGroup.FEM.Tests.Elements
{
	/// <summary>
	/// Tests 8-noded quadrilateral instances of <see cref="ContinuumElement2D"/> against Abaqus.
	/// Authors: Serafeim Bakalakos
	/// </summary>
	public class Quad8
	{
		private static double thickness = 1.0;

		private static readonly ElasticMaterial2D material0 = new ElasticMaterial2D(StressState2D.PlaneStress)
		{
			YoungModulus = 2.1e5,
			PoissonRatio = 0.3
		};

		private static readonly DynamicMaterial dynamicMaterial = new DynamicMaterial(78.5, 0, 0);

		/// <summary>
		/// Random shape, not too distorted.
		/// </summary>
		private static readonly IReadOnlyList<Node> nodeSet0 = new Node[]
		{
			new Node( id: 0, x: 0.7, y:  2.0 ),
			new Node( id: 1, x: 0.2, y:  0.3 ),
			new Node( id: 2, x: 2.0, y:  0.9 ),
			new Node( id: 3, x: 3.0, y:  2.7 ),

			new Node( id: 4, x: 0.7, y:  1.1 ),
			new Node( id: 5, x: 1.3, y:  0.1 ),
			new Node( id: 6, x: 2.1, y:  1.9 ),
			new Node( id: 7, x: 1.8, y:  2.5 ),
		};

		/// <summary>
		/// To reproduce the reference solution, run the Abaqus job quad8_test0.inp and look at 
		/// QUAD8_TEST0_MASS_MATRICES.mtx.
		/// </summary>
		[Fact]
		private static void TestConsistentMass0()
		{
			IQuadrature2D quadratureForMass = GaussLegendre2D.GetQuadratureWithOrder(3, 3);

			var materialsAtGaussPoints = new List<ElasticMaterial2D>();
			foreach (GaussPoint gaussPoint in quadratureForMass.IntegrationPoints)
			{
				materialsAtGaussPoints.Add(material0.Clone());
			}
			var quad8 = new ContinuumElement2D(thickness, nodeSet0, InterpolationQuad8.UniqueInstance,
				GaussLegendre2D.GetQuadratureWithOrder(3, 3), quadratureForMass,
				ExtrapolationGaussLegendre3x3.UniqueInstance,
				materialsAtGaussPoints, dynamicMaterial);

			IMatrix M = quad8.BuildConsistentMassMatrix();
			Matrix expectedM = Matrix.CreateFromArray(new double[,]
			{
				{ 9.115245555556, 0.000000000000, 1.881383333333, 0.000000000000, 3.914882222222, 0.000000000000, 2.417800000000, 0.000000000000, -7.376906666667, 0.000000000000, -11.056637777778, 0.000000000000, -9.104604444445, 0.000000000000, -7.444940000000, 0.000000000000      },
				{ 0.000000000000, 9.115245555556, 0.000000000000, 1.881383333333, 0.000000000000, 3.914882222222, 0.000000000000, 2.417800000000, 0.000000000000, -7.376906666667, 0.000000000000, -11.056637777778, 0.000000000000, -9.104604444445, 0.000000000000, -7.444940000000      },
				{ 1.881383333333, 0.000000000000, 9.727545555556, 0.000000000000, 2.239866666667, 0.000000000000, 3.841615555556, 0.000000000000, -6.727973333334, 0.000000000000, -6.701806666667, 0.000000000000, -9.031337777778, 0.000000000000, -11.077571111111, 0.000000000000      },
				{ 0.000000000000, 1.881383333333, 0.000000000000, 9.727545555556, 0.000000000000, 2.239866666667, 0.000000000000, 3.841615555556, 0.000000000000, -6.727973333334, 0.000000000000, -6.701806666667, 0.000000000000, -9.031337777778, 0.000000000000, -11.077571111111      },
				{ 3.914882222222, 0.000000000000, 2.239866666667, 0.000000000000, 7.681312222222, 0.000000000000, 2.776283333333, 0.000000000000, -11.082804444445, 0.000000000000, -8.784673333334, 0.000000000000, -6.832640000000, 0.000000000000, -11.150837777778, 0.000000000000     },
				{ 0.000000000000, 3.914882222222, 0.000000000000, 2.239866666667, 0.000000000000, 7.681312222222, 0.000000000000, 2.776283333333, 0.000000000000, -11.082804444445, 0.000000000000, -8.784673333334, 0.000000000000, -6.832640000000, 0.000000000000, -11.150837777778     },
				{ 2.417800000000, 0.000000000000, 3.841615555556, 0.000000000000, 2.776283333333, 0.000000000000, 7.581878888889, 0.000000000000, -11.009537777778, 0.000000000000, -10.983371111111, 0.000000000000, -6.895440000000, 0.000000000000, -8.941673333334, 0.000000000000     },
				{ 0.000000000000, 2.417800000000, 0.000000000000, 3.841615555556, 0.000000000000, 2.776283333333, 0.000000000000, 7.581878888889, 0.000000000000, -11.009537777778, 0.000000000000, -10.983371111111, 0.000000000000, -6.895440000000, 0.000000000000, -8.941673333334     },
				{ -7.376906666667, 0.000000000000, -6.727973333334, 0.000000000000, -11.082804444445, 0.000000000000, -11.009537777778, 0.000000000000, 45.023413333334, 0.000000000000, 28.692622222223, 0.000000000000, 20.114142222223, 0.000000000000, 28.734488888889, 0.000000000000 },
				{ 0.000000000000, -7.376906666667, 0.000000000000, -6.727973333334, 0.000000000000, -11.082804444445, 0.000000000000, -11.009537777778, 0.000000000000, 45.023413333334, 0.000000000000, 28.692622222223, 0.000000000000, 20.114142222223, 0.000000000000, 28.734488888889 },
				{ -11.056637777778, 0.000000000000, -6.701806666667, 0.000000000000, -8.784673333334, 0.000000000000, -10.983371111111, 0.000000000000, 28.692622222223, 0.000000000000, 48.215746666667, 0.000000000000, 24.589688888889, 0.000000000000, 22.134208888889, 0.000000000000 },
				{ 0.000000000000, -11.056637777778, 0.000000000000, -6.701806666667, 0.000000000000, -8.784673333334, 0.000000000000, -10.983371111111, 0.000000000000, 28.692622222223, 0.000000000000, 48.215746666667, 0.000000000000, 24.589688888889, 0.000000000000, 22.134208888889 },
				{ -9.104604444445, 0.000000000000, -9.031337777778, 0.000000000000, -6.832640000000, 0.000000000000, -6.895440000000, 0.000000000000, 20.114142222223, 0.000000000000, 24.589688888889, 0.000000000000, 28.821013333334, 0.000000000000, 24.924622222223, 0.000000000000   },
				{ 0.000000000000, -9.104604444445, 0.000000000000, -9.031337777778, 0.000000000000, -6.832640000000, 0.000000000000, -6.895440000000, 0.000000000000, 20.114142222223, 0.000000000000, 24.589688888889, 0.000000000000, 28.821013333334, 0.000000000000, 24.924622222223   },
				{ -7.444940000000, 0.000000000000, -11.077571111111, 0.000000000000, -11.150837777778, 0.000000000000, -8.941673333334, 0.000000000000, 28.734488888889, 0.000000000000, 22.134208888889, 0.000000000000, 24.924622222223, 0.000000000000, 49.869480000000, 0.000000000000 },
				{ 0.000000000000, -7.444940000000, 0.000000000000, -11.077571111111, 0.000000000000, -11.150837777778, 0.000000000000, -8.941673333334, 0.000000000000, 28.734488888889, 0.000000000000, 22.134208888889, 0.000000000000, 24.924622222223, 0.000000000000, 49.869480000000 }
			}); // from Abaqus
			Assert.True(M.Equals(expectedM, 1e-10));
		}

		/// <summary>
		/// To reproduce the reference solution, run the Abaqus job quad8_test0.inp and look at 
		/// QUAD8_TEST0_STIFFNESS_MATRICES.mtx.
		/// </summary>
		[Fact]
		private static void TestStiffness0()
		{
			var factory = new ContinuumElement2DFactory(thickness, material0, dynamicMaterial);
			ContinuumElement2D quad8 = factory.CreateElement(CellType.Quad8, nodeSet0);
			IMatrix K = quad8.BuildStiffnessMatrix();
			double[,] expectedK = new double[,]
			{
				{ 375957.1179055000, -130889.8953535400, 150398.7276607000, -30577.4536395990, 256263.9470137800, -74238.0130438120, 120262.5252654800, -29264.7336935000, -201002.5528104600, 86903.8674540020, -257013.8069723100, 62503.0900498910, -46878.8680758500, 18155.7769717120, -397987.0899868400, 97407.3612548510    },
				{ -130889.8953535400, 264620.9765739100, -32500.5305626760, 108395.9237548000, -74238.0130438120, 148536.2310068900, -27341.6567704230, 66761.0027359010, 94596.1751463100, -240343.8078467800, 62503.0900498910, -121948.9418401700, 18155.7769717120, -45778.2662465080, 89715.0535625430, -180243.1181380300     },
				{ 150398.7276607000, -32500.5305626760, 227070.8522891000, 24999.1791433860, 228075.7272750800, -35916.2249916640, 84715.8706078410, 5691.9746958972, -125567.3958918900, -7749.8459233937, -292871.0398821500, 29891.4228662980, -48734.4768523870, -12619.4433884890, -223088.2652062900, 28203.4681606430        },
				{ -30577.4536395990, 108395.9237548000, 24999.1791433860, 135717.9797677600, -37839.3019147410, 115695.2678105500, 5691.9746958972, 47395.7262976820, -15442.1536157010, -174789.3743509800, 37583.7305586050, -85742.6925983560, -12619.4433884890, -59359.2032811550, 28203.4681606430, -87313.6274002970         },
				{ 256263.9470137800, -74238.0130438120, 228075.7272750800, -37839.3019147410, 758698.8330168600, -151201.1982884900, 164115.4999410500, -25847.0928360870, -256843.4834062300, 67264.1450855650, -543079.2639940800, 68308.1044726010, -237729.8010194900, 72706.8336500980, -369501.4588269600, 80846.5228748690   },
				{ -74238.0130438120, 148536.2310068900, -35916.2249916640, 115695.2678105500, -151201.1982884900, 368478.4546134000, -27770.1697591640, 85552.6038287820, 67264.1450855650, -166187.6357417500, 60615.7967802930, -223855.9737213700, 80399.1413424050, -193251.7815799700, 80846.5228748690, -134967.1662165300    },
				{ 120262.5252654800, -27341.6567704230, 84715.8706078410, 5691.9746958972, 164115.4999410500, -27770.1697591640, 120174.7405891600, 2419.9582040241, -64973.3681564819, 9886.9381472210, -156296.6219667900, 10235.6574488430, -29392.1129538440, -12237.5732893820, -238606.5333264100, 39114.8713229830           },
				{ -29264.7336935000, 66761.0027359010, 5691.9746958972, 47395.7262976820, -25847.0928360870, 85552.6038287820, 2419.9582040241, 63320.7521890160, 9886.9381472210, -67395.8289291320, 10235.6574488430, -47867.5471489640, -19929.8809816900, -81314.4790110740, 46807.1790152900, -66452.2299622090                },
				{ -201002.5528104600, 94596.1751463100, -125567.3958918900, -15442.1536157010, -256843.4834062300, 67264.1450855650, -64973.3681564819, 9886.9381472210, 429358.6456627000, -85681.3574299900, 91327.9338391870, 31186.8509038820, -120846.2914647200, 7660.6709468983, 248546.5122279100, -109471.2691841800       },
				{ 86903.8674540020, -240343.8078467800, -7749.8459233937, -174789.3743509800, 67264.1450855650, -166187.6357417500, 9886.9381472210, -67395.8289291320, -85681.3574299900, 501789.8020087000, 31186.8509038820, -24500.7443298540, 7660.6709468983, 60837.6092548520, -109471.2691841800, 110589.9799349400         },
				{ -257013.8069723100, 62503.0900498910, -292871.0398821500, 37583.7305586050, -543079.2639940800, 60615.7967802930, -156296.6219667900, 10235.6574488430, 91327.9338391870, 31186.8509038820, 618909.3895160900, -107045.7505943100, 156953.5398894400, -38863.3219066200, 382069.8695706200, -56216.0532405820     },
				{ 62503.0900498910, -121948.9418401700, 29891.4228662980, -85742.6925983560, 68308.1044726010, -223855.9737213700, 10235.6574488430, -47867.5471489640, 31186.8509038820, -24500.7443298540, -107045.7505943100, 363810.1766979700, -38863.3219066200, 60473.5634510330, -56216.0532405820, 79632.1594897160        },
				{ -46878.8680758500, 18155.7769717120, -48734.4768523870, -12619.4433884890, -237729.8010194900, 80399.1413424050, -29392.1129538440, -19929.8809816900, -120846.2914647200, 7660.6709468983, 156953.5398894400, -38863.3219066200, 585811.6805262000, -190782.6753367700, -259183.6700493400, 155979.7323525500    },
				{ 18155.7769717120, -45778.2662465080, -12619.4433884890, -59359.2032811550, 72706.8336500980, -193251.7815799700, -12237.5732893820, -81314.4790110740, 7660.6709468983, 60837.6092548520, -38863.3219066200, 60473.5634510330, -190782.6753367700, 616593.8061240300, 155979.7323525500, -358201.2487112100       },
				{ -397987.0899868400, 89715.0535625430, -223088.2652062900, 28203.4681606430, -369501.4588269600, 80846.5228748690, -238606.5333264100, 46807.1790152900, 248546.5122279100, -109471.2691841800, 382069.8695706200, -56216.0532405820, -259183.6700493400, 155979.7323525500, 857750.6355973100, -235864.6335411300 },
				{ 97407.3612548510, -180243.1181380300, 28203.4681606430, -87313.6274002970, 80846.5228748690, -134967.1662165300, 39114.8713229830, -66452.2299622090, -109471.2691841800, 110589.9799349400, -56216.0532405820, 79632.1594897160, 155979.7323525500, -358201.2487112100, -235864.6335411300, 636955.2510036100    }
			}; // from Abaqus
			Assert.True(K.Equals(Matrix.CreateFromArray(expectedK), 1e-10));
		}

		/// <summary>
		/// To reproduce the reference solution, run the Abaqus job quad8_test0.inp and look at tri8_test0.dat.
		/// nodes).
		/// </summary>
		[Fact]
		public static void TestStrainsStresses0()
		{
			var factory = new ContinuumElement2DFactory(thickness, material0, null);
			ContinuumElement2D quad8 = factory.CreateElement(CellType.Quad8, nodeSet0);

			// Abaqus results
			double[] displacements =
			{
				0.0, 0.0,                   // Node 1
                0.0, 0.0,                   // Node 2
                -6.1872E-03, -1.3677E-02,   // Node 3
                2.1193E-02, -6.0952E-02,    // Node 4
                -6.4016E-03, -2.9598E-03,   // Node 5
                -9.8052E-03, -7.6739E-03,   // Node 6
                2.0961E-03, -2.4904E-02,    // Node 7
                1.2842E-02, -2.0029E-02     // Node 8
            };

			double[][] expectedStrainsAtGPs =
			{
				new double[] {  5.0483E-03,  9.3636E-04,  1.3684E-03 },  // Gauss point 1
                new double[] {  4.4925E-04,  6.0911E-04, -1.1381E-03 },  // Gauss point 2
                new double[] { -1.0021E-02,  3.7130E-04, -6.7670E-03 },  // Gauss point 3
                new double[] {  3.0367E-03, -4.2664E-03, -6.6618E-03 },  // Gauss point 4
                new double[] {  7.6921E-04, -2.1368E-03, -5.1775E-03 },  // Gauss point 5
                new double[] { -2.5482E-03,  1.8783E-03, -6.8114E-03 },  // Gauss point 6
                new double[] {  1.2133E-03, -2.1843E-05, -1.4123E-02 },  // Gauss point 7
                new double[] { -4.6295E-03, -1.0028E-02, -4.5072E-03 },  // Gauss point 8
                new double[] {  4.6685E-03, -2.5173E-03, -3.6237E-03 }   // Gauss point 9
            };
			double[][] expectedStressesAtGPs =
			{
				new double[] {  1230.00,   565.60,    110.50 },  // Gauss point 1
                new double[] {   145.80,   171.70,    -91.92 },  // Gauss point 2
                new double[] { -2287.00,  -608.10,   -546.60 },  // Gauss point 3
                new double[] {   405.40,  -774.30,   -538.10 },  // Gauss point 4
                new double[] {    29.58,  -439.80,   -418.20 },  // Gauss point 5
                new double[] {  -458.00,   257.00,   -550.20 },  // Gauss point 6
                new double[] {   278.50,    78.95,  -1141.00 },  // Gauss point 7
                new double[] { -1763.00, -2635.00,   -364.00 },  // Gauss point 8
                new double[] {   903.10,  -257.70,   -292.70 }  // Gauss point 9
            };

			// The order of the nodes is the same (at least in this job)
			double[][] expectedStrainsAtNodes =
			{
				new double[] {  6.5625E-03,  5.7325E-03,  3.8681E-03 },  // Node 1
                new double[] { -1.6562E-02, -1.2399E-03, -8.9890E-03 },  // Node 2
                new double[] {  1.4923E-02,  2.1676E-03, -3.8448E-03 },  // Node 3
                new double[] {  7.1739E-03,  1.2357E-02, -2.2045E-02 },  // Node 4
                new double[] { -7.1803E-04,  4.4174E-04,  9.2207E-04 },  // Node 5
                new double[] { -3.7108E-03,  3.4008E-03, -7.8726E-03 },  // Node 6
                new double[] { -7.2747E-03, -1.3290E-02, -3.4274E-03 },  // Node 7
                new double[] {  3.4993E-03, -4.5319E-03, -7.6795E-03 }  // Node 8
            };
			double[][] expectedStressesAtNodes =
			{
				new double[] {  1911.0,  1777.00,   312.40 },  // Node 1
                new double[] { -3908.0, -1433.00,  -726.00 },  // Node 2
                new double[] {  3594.0,  1533.00,  -310.50 },  // Node 3
                new double[] {  2511.0,  3348.00, -1781.00 },  // Node 4
                new double[] {  -135.1,    52.23,    74.47 },  // Node 5
                new double[] {  -620.9,   527.90,  -635.90 },  // Node 6
                new double[] { -2599.0, -3571.00,  -276.80 },  // Node 7
                new double[] {   493.8, - 803.60,  -620.30 }   // Node 8
            };

			(IReadOnlyList<double[]> strainsAtGPs, IReadOnlyList<double[]> stressesAtGPs) =
				quad8.UpdateStrainsStressesAtGaussPoints(displacements);
			IReadOnlyList<double[]> strainsAtNodes =
				quad8.GaussPointExtrapolation.ExtrapolateTensorFromGaussPointsToNodes(strainsAtGPs, quad8.Interpolation);
			IReadOnlyList<double[]> stressesAtNodes =
				quad8.GaussPointExtrapolation.ExtrapolateTensorFromGaussPointsToNodes(stressesAtGPs, quad8.Interpolation);

			Assert.True(Utilities.AreTensorsEqual(expectedStrainsAtGPs, strainsAtGPs, 1e-3));
			Assert.True(Utilities.AreTensorsEqual(expectedStressesAtGPs, stressesAtGPs, 0.5e-2));
			Assert.True(Utilities.AreTensorsEqual(expectedStrainsAtNodes, strainsAtNodes, 1e-3));
			Assert.True(Utilities.AreTensorsEqual(expectedStressesAtNodes, stressesAtNodes, 0.5e-2));
		}
	}
}
