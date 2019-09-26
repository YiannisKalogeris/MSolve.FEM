﻿using System.Collections.Generic;
using MGroup.FEM.Elements;
using MGroup.FEM.Entities;
using MGroup.LinearAlgebra.Matrices;
using MGroup.Materials;
using MGroup.MSolve.Discretization.Mesh;
using Xunit;

namespace MGroup.FEM.Tests.Elements
{
	/// <summary>
	/// Tests 8-noded hexahedral instances of <see cref="Hexa8"/>
	/// </summary>
	public class Hexa8
	{
		private static readonly ElasticMaterial3D Material0 = new ElasticMaterial3D
		{
			YoungModulus = 1000,
			PoissonRatio = 0.3
		};

		private static readonly ElasticMaterial3D Material1 = new ElasticMaterial3D
		{
			YoungModulus = 210000,
			PoissonRatio = 0.3
		};

		private static readonly DynamicMaterial DynamicMaterial0 = new DynamicMaterial(1, 0, 0);
		private static readonly DynamicMaterial DynamicMaterial1 = new DynamicMaterial(78.5, 0, 0);

		// The nodes in the order of the Abaqus Hexa8 element
		private static readonly IReadOnlyList<Node> NodeSet0 = new Node[]
		{
			new Node( id: 0, x: 0, y:  5, z: 5 ),
			new Node( id: 1, x: 0, y:  0, z: 5 ),
			new Node( id: 2, x: 0, y:  0, z: 0 ),
			new Node( id: 3, x: 0, y:  5, z: 0 ),

			new Node( id: 4, x: 5, y:  5, z: 5 ),
			new Node( id: 5, x: 5, y:  0, z: 5 ),
			new Node( id: 6, x: 5, y:  0, z: 0 ),
			new Node( id: 7, x: 5, y:  5, z: 0 ),
		};

		// The nodes in the order of the Abaqus Hexa8 element
		private static readonly IReadOnlyList<Node> NodeSet1 = new Node[]
		{
			new Node( id: 0, x: -1, y:   1, z: 2 ),
			new Node( id: 1, x: -1, y:  -1, z: 2 ),
			new Node( id: 2, x: -1, y:  -1, z: 0 ),
			new Node( id: 3, x: -1, y:   1, z: 0 ),

			new Node( id: 4, x:  1, y:   1, z: 2 ),
			new Node( id: 5, x:  1, y:  -1, z: 2 ),
			new Node( id: 6, x:  1, y:  -1, z: 0 ),
			new Node( id: 7, x:  1, y:   1, z: 0 ),
		};

		internal static IReadOnlyList<Node> AbaqusToMSolveHexa8Nodes(IReadOnlyList<Node> abaqusOrder)
		{
			var msolveFromAbaqusMap = new int[] { 2, 6, 7, 3, 1, 5, 4, 0 };
			var msolveOrder = new Node[8];
			for (int i = 0; i < 8; ++i) msolveOrder[i] = abaqusOrder[msolveFromAbaqusMap[i]];
			return msolveOrder;
		}

		internal static Matrix AbaqusToMSolveHexa8StructuralMatrix(IMatrixView abaqusMatrix)
		{
			var msolveFromAbaqusNodeMap = new int[] { 2, 6, 7, 3, 1, 5, 4, 0 };
			var msolveFromAbaqusDofMap = new int[24];
			for (int node = 0; node < 8; ++node)
			{
				msolveFromAbaqusDofMap[3 * node + 0] = 3 * msolveFromAbaqusNodeMap[node] + 0;
				msolveFromAbaqusDofMap[3 * node + 1] = 3 * msolveFromAbaqusNodeMap[node] + 1;
				msolveFromAbaqusDofMap[3 * node + 2] = 3 * msolveFromAbaqusNodeMap[node] + 2;
			}

			Matrix msolveMatrix = abaqusMatrix.CopyToFullMatrix();
			return msolveMatrix.Reorder(msolveFromAbaqusDofMap, false);
		}

		//[Fact]
		private static void TestStiffnessMatrix0()
		{
			var factory = new ContinuumElement3DFactory(Material0, DynamicMaterial0);
			var hexa8 = factory.CreateElement(CellType.Hexa8, AbaqusToMSolveHexa8Nodes(NodeSet0));
			var computedK = hexa8.BuildStiffnessMatrix();

			var expectedK = AbaqusToMSolveHexa8StructuralMatrix(Matrix.CreateFromArray(new double[24, 24]
			{
				{ 972.667378917380  ,  -313.835470085470 ,  -313.835470085470 ,  296.029202279200  ,  -6.67735042735140  , -287.126068376070 ,  296.029202279200  ,  -287.126068376070 ,  -6.67735042735140 ,  117.966524216520  ,  126.869658119660  ,  126.869658119660  ,  -331.641737891740 ,  6.67735042735140 ,   6.67735042735140  ,  -456.285612535610 ,  313.835470085470  ,  -126.869658119660 ,  -456.285612535610 ,  -126.869658119660 ,  313.835470085470  ,  -438.479344729340 ,  287.126068376070   , 287.126068376070},
				{ -313.835470085470 ,  972.667378917380 ,   313.835470085470  ,  6.67735042735140  ,  -331.641737891740  , -6.67735042735140 ,  -287.126068376070 ,  296.029202279200  ,  6.67735042735140 ,   -126.869658119660 ,  -456.285612535610 ,  -313.835470085470 ,  -6.67735042735140 ,  296.029202279200 ,   287.126068376070  ,  313.835470085470  ,  -456.285612535610 ,  126.869658119660  ,  126.869658119660  ,  117.966524216520  ,  -126.869658119660 ,  287.126068376070  ,  -438.479344729340  , -287.126068376070},
				{ -313.835470085470 ,  313.835470085470 ,   972.667378917380  ,  -287.126068376070 ,  6.67735042735140   , 296.029202279200  ,  6.67735042735140  ,  -6.67735042735140 ,  -331.641737891740,   -126.869658119660 ,  -313.835470085470 ,  -456.285612535610 ,  -6.67735042735140 ,  287.126068376070 ,   296.029202279200  ,  126.869658119660  ,  -126.869658119660 ,  117.966524216520  ,  313.835470085470  ,  126.869658119660  ,  -456.285612535610 ,  287.126068376070  ,  -287.126068376070  , -438.479344729340},
				{ 296.029202279200  ,  6.67735042735140  ,  -287.126068376070 ,  972.667378917380  ,  313.835470085470  ,  -313.835470085470 ,  117.966524216520  ,  -126.869658119660 ,  126.869658119660 ,   296.029202279200  ,  287.126068376070  ,  -6.67735042735140 ,  -456.285612535610 ,  -313.835470085470,   -126.869658119660 ,  -331.641737891740 ,  -6.67735042735140 ,  6.67735042735140  ,  -438.479344729340 ,  -287.126068376070 ,  287.126068376070  ,  -456.285612535610 ,  126.869658119660   , 313.835470085470},
				{ -6.67735042735140 ,  -331.641737891740 ,  6.67735042735140  ,  313.835470085470  ,  972.667378917380  ,  -313.835470085470 ,  126.869658119660  ,  -456.285612535610 ,  313.835470085470 ,   287.126068376070  ,  296.029202279200  ,  -6.67735042735140,   -313.835470085470 ,  -456.285612535610,   -126.869658119660 ,  6.67735042735130  ,  296.029202279200 ,   -287.126068376070 ,  -287.126068376070 ,  -438.479344729340 ,  287.126068376070  ,  -126.869658119660 ,  117.966524216520   , 126.869658119660},
				{ -287.126068376070 ,  -6.67735042735140 ,  296.029202279200  ,  -313.835470085470 ,  -313.835470085470 ,  972.667378917380  ,  -126.869658119660 ,  313.835470085470  ,  -456.285612535610,   6.67735042735140  ,  6.67735042735140  ,  -331.641737891740,   126.869658119660  ,  126.869658119660 ,   117.966524216520  ,  -6.67735042735140 ,  -287.126068376070,   296.029202279200  ,  287.126068376070  ,  287.126068376070  ,  -438.479344729340 ,  313.835470085470  ,  -126.869658119660  , -456.285612535610},
				{ 296.029202279200  ,  -287.126068376070 ,  6.67735042735140  ,  117.966524216520  ,  126.869658119660  ,  -126.869658119660 ,  972.667378917380  ,  -313.835470085470 ,  313.835470085470,    296.029202279200  ,  -6.67735042735140 ,  287.126068376070 ,   -456.285612535610 ,  -126.869658119660,   -313.835470085470 ,  -438.479344729340 ,  287.126068376070 ,   -287.126068376070 ,  -331.641737891740 ,  6.67735042735140  ,  -6.67735042735140 ,  -456.285612535610 ,  313.835470085470   , 126.869658119660},
				{ -287.126068376070 ,  296.029202279200   , -6.67735042735140 ,  -126.869658119660 ,  -456.285612535610 ,  313.835470085470  ,  -313.835470085470 ,  972.667378917380  ,  -313.835470085470 ,  6.67735042735140  ,  -331.641737891740 ,  6.67735042735130  ,  126.869658119660  ,  117.966524216520 ,   126.869658119660  ,  287.126068376070  ,  -438.479344729340 ,  287.126068376070  ,  -6.67735042735140 ,  296.029202279200  ,  -287.126068376070 ,  313.835470085470  ,  -456.285612535610  , -126.869658119660},
				{ -6.67735042735140 ,  6.67735042735140  ,  -331.641737891740 ,  126.869658119660  ,  313.835470085470  ,  -456.285612535610 ,  313.835470085470  ,  -313.835470085470 ,  972.667378917380  ,  287.126068376070  ,  -6.67735042735140 ,  296.029202279200  ,  -313.835470085470 ,  -126.869658119660,   -456.285612535610 ,  -287.126068376070 ,  287.126068376070  ,  -438.479344729340 ,  6.67735042735140  ,  -287.126068376070 ,  296.029202279200  ,  -126.869658119660 ,  126.869658119660   , 117.966524216520},
				{ 117.966524216520  ,  -126.869658119660 ,  -126.869658119660 ,  296.029202279200  ,  287.126068376070  ,  6.67735042735140  ,  296.029202279200  ,  6.67735042735140  ,  287.126068376070  ,  972.667378917380  ,  313.835470085470  ,  313.835470085470  ,  -438.479344729340 ,  -287.126068376070,   -287.126068376070 ,  -456.285612535610 ,  126.869658119660  ,  -313.835470085470 ,  -456.285612535610 ,  -313.835470085470 ,  126.869658119660  ,  -331.641737891740 ,  -6.67735042735140  , -6.67735042735140},
				{ 126.869658119660   , -456.285612535610  , -313.835470085470 ,  287.126068376070  ,  296.029202279200  ,  6.67735042735140  ,  -6.67735042735140 ,  -331.641737891740 ,  -6.67735042735140 ,  313.835470085470  ,  972.667378917380  ,  313.835470085470  ,  -287.126068376070 ,  -438.479344729340,   -287.126068376070 ,  -126.869658119660 ,  117.966524216520  ,  -126.869658119660 ,  -313.835470085470 ,  -456.285612535610 ,  126.869658119660  ,  6.67735042735140  ,  296.029202279200   , 287.126068376070},
				{ 126.869658119660  ,  -313.835470085470 ,  -456.285612535610 ,  -6.67735042735140 ,  -6.67735042735140 ,  -331.641737891740 ,  287.126068376070  ,  6.67735042735130  ,  296.029202279200  ,  313.835470085470  ,  313.835470085470  ,  972.667378917380  ,  -287.126068376070 ,  -287.126068376070,   -438.479344729340 ,  -313.835470085470 ,  126.869658119660  ,  -456.285612535610 ,  -126.869658119660  , -126.869658119660 ,  117.966524216520  ,  6.67735042735140  ,  287.126068376070   , 296.029202279200},
				{ -331.641737891740 ,  -6.67735042735140  , -6.67735042735140 ,  -456.285612535610 ,  -313.835470085470 ,  126.869658119660  ,  -456.285612535610 ,  126.869658119660  ,  -313.835470085470 ,  -438.479344729340 ,  -287.126068376070 ,  -287.126068376070 ,  972.667378917380  ,  313.835470085470 ,   313.835470085470  ,  296.029202279200  ,  6.67735042735130  ,  287.126068376070  ,  296.029202279200  ,  287.126068376070  ,  6.67735042735140  ,  117.966524216520  ,  -126.869658119660  , -126.869658119660},
				{ 6.67735042735140  ,  296.029202279200  ,  287.126068376070  ,  -313.835470085470 ,  -456.285612535610 ,  126.869658119660  ,  -126.869658119660 ,  117.966524216520  ,  -126.869658119660 ,  -287.126068376070 ,  -438.479344729340 ,  -287.126068376070 ,  313.835470085470  ,  972.667378917380 ,   313.835470085470  ,  -6.67735042735140 ,  -331.641737891740 ,  -6.67735042735140 ,  287.126068376070  ,  296.029202279200  ,  6.67735042735140  ,  126.869658119660  ,  -456.285612535610  , -313.835470085470},
				{ 6.67735042735140  ,  287.126068376070  ,  296.029202279200  ,  -126.869658119660 ,  -126.869658119660 ,  117.966524216520  ,  -313.835470085470 ,  126.869658119660  ,  -456.285612535610 ,  -287.126068376070 ,  -287.126068376070 ,  -438.479344729340 ,  313.835470085470  ,  313.835470085470 ,   972.667378917380  ,  287.126068376070  ,  6.67735042735140  ,  296.029202279200  ,  -6.67735042735140 ,  -6.67735042735140 ,  -331.641737891740 ,  126.869658119660  ,  -313.835470085470  , -456.285612535610},
				{ -456.285612535610 ,  313.835470085470  ,  126.869658119660  ,  -331.641737891740 ,  6.67735042735130  ,  -6.67735042735140 ,  -438.479344729340 ,  287.126068376070  ,  -287.126068376070 ,  -456.285612535610 ,  -126.869658119660 ,  -313.835470085470 ,  296.029202279200  ,  -6.67735042735140,   287.126068376070  ,  972.667378917380  ,  -313.835470085470 ,  313.835470085470  ,  117.966524216520  ,  126.869658119660  ,  -126.869658119660 ,  296.029202279200  ,  -287.126068376070  , 6.67735042735140},
				{ 313.835470085470  ,  -456.285612535610 ,  -126.869658119660 ,  -6.67735042735140 ,  296.029202279200  ,  -287.126068376070 ,  287.126068376070  ,  -438.479344729340 ,  287.126068376070  ,  126.869658119660  ,  117.966524216520  ,  126.869658119660  ,  6.67735042735130  ,  -331.641737891740,   6.67735042735140  ,  -313.835470085470 ,  972.667378917380  ,  -313.835470085470 ,  -126.869658119660 ,  -456.285612535610 ,  313.835470085470  ,  -287.126068376070 ,  296.029202279200   , -6.67735042735140},
				{ -126.869658119660 ,  126.869658119660  ,  117.966524216520  ,  6.67735042735140  ,  -287.126068376070 ,  296.029202279200  ,  -287.126068376070 ,  287.126068376070  ,  -438.479344729340 ,  -313.835470085470 ,  -126.869658119660 ,  -456.285612535610 ,  287.126068376070  ,  -6.67735042735140,   296.029202279200  ,  313.835470085470  ,  -313.835470085470 ,  972.667378917380  ,  126.869658119660  ,  313.835470085470  ,  -456.285612535610 ,  -6.67735042735130 ,  6.67735042735140   , -331.641737891740},
				{ -456.285612535610  , 126.869658119660  ,  313.835470085470  ,  -438.479344729340 ,  -287.126068376070 ,  287.126068376070  ,  -331.641737891740 ,  -6.67735042735140 ,  6.67735042735140  ,  -456.285612535610 ,  -313.835470085470 ,  -126.869658119660 ,  296.029202279200  ,  287.126068376070 ,   -6.67735042735140 ,  117.966524216520  ,  -126.869658119660 ,  126.869658119660  ,  972.667378917380  ,  313.835470085470  ,  -313.835470085470 ,  296.029202279200  ,  6.67735042735140   , -287.126068376070},
				{ -126.869658119660 ,  117.966524216520  ,  126.869658119660  ,  -287.126068376070 ,  -438.479344729340 ,  287.126068376070  ,  6.67735042735140  ,  296.029202279200  ,  -287.126068376070 ,  -313.835470085470 ,  -456.285612535610 ,  -126.869658119660 ,  287.126068376070  ,  296.029202279200 ,   -6.67735042735140 ,  126.869658119660  ,  -456.285612535610 ,  313.835470085470  ,  313.835470085470  ,  972.667378917380  ,  -313.835470085470 ,  -6.67735042735140 ,  -331.641737891740  , 6.67735042735140},
				{ 313.835470085470  ,  -126.869658119660 ,  -456.285612535610 ,  287.126068376070  ,  287.126068376070  ,  -438.479344729340,   -6.67735042735140 ,  -287.126068376070 ,  296.029202279200  ,  126.869658119660  ,  126.869658119660  ,  117.966524216520  ,  6.67735042735140  ,  6.67735042735140 ,   -331.641737891740 ,  -126.869658119660 ,  313.835470085470  ,  -456.285612535610 ,  -313.835470085470 ,  -313.835470085470 ,  972.667378917380  ,  -287.126068376070 ,  -6.67735042735140  , 296.029202279200},
				{ -438.479344729340 ,  287.126068376070  ,  287.126068376070  ,  -456.285612535610 ,  -126.869658119660 ,  313.835470085470  ,  -456.285612535610 ,  313.835470085470  ,  -126.869658119660 ,  -331.641737891740 ,  6.67735042735140  ,  6.67735042735140  ,  117.966524216520  ,  126.869658119660 ,   126.869658119660  ,  296.029202279200  ,  -287.126068376070 ,  -6.67735042735130 ,  296.029202279200  ,  -6.67735042735140 ,  -287.126068376070 ,  972.667378917380  ,  -313.835470085470  , -313.835470085470},
				{ 287.126068376070  ,  -438.479344729340 ,  -287.126068376070 ,  126.869658119660  ,  117.966524216520  ,  -126.869658119660 ,  313.835470085470  ,  -456.285612535610 ,  126.869658119660  ,  -6.67735042735140 ,  296.029202279200  ,  287.126068376070  ,  -126.869658119660 ,  -456.285612535610,   -313.835470085470 ,  -287.126068376070 ,  296.029202279200  ,  6.67735042735140  ,  6.67735042735140  ,  -331.641737891740 ,  -6.67735042735140 ,  -313.835470085470 ,  972.667378917380   , 313.835470085470},
				{ 287.126068376070   , -287.126068376070  , -438.479344729340 ,  313.835470085470  ,  126.869658119660  ,  -456.285612535610 ,  126.869658119660  ,  -126.869658119660 ,  117.966524216520  ,  -6.67735042735140 ,  287.126068376070  ,  296.029202279200  ,  -126.869658119660 ,  -313.835470085470,   -456.285612535610 ,  6.67735042735140  ,  -6.67735042735140 ,  -331.641737891740 ,  -287.126068376070 ,  6.67735042735140  ,  296.029202279200  ,  -313.835470085470 ,  313.835470085470   , 972.667378917380},
			}));

			Assert.True(computedK.Equals(expectedK, 1e-8));
		}

		/// <summary>
		/// To reproduce the reference solution, run the Abaqus job hexa8_cube.inp and look at 
		/// HEXA8_CUBE_STIFFNESS_MATRICES.mtx
		/// </summary>
		//[Fact]
		private static void TestStiffnessMatrix1()
		{
			var factory = new ContinuumElement3DFactory(Material1, DynamicMaterial1);
			var hexa8 = factory.CreateElement(CellType.Hexa8, AbaqusToMSolveHexa8Nodes(NodeSet1));
			var computedK = hexa8.BuildStiffnessMatrix();

			//string output = @"C:\Users\Serafeim\Desktop\hexa8.txt";
			//var writer = new LinearAlgebra.Output.FullMatrixWriter();
			////writer.ArrayFormat = LinearAlgebra.Output.Formatting.Array2DFormat.CSharpArray;
			//writer.WriteToFile(computedK, output);

			var expectedK = AbaqusToMSolveHexa8StructuralMatrix(Matrix.CreateFromArray(new double[24, 24]
			{
				{ 81704.05983, -26362.17949, -26362.17949, 24866.45299, -560.8974359, -24118.58974, 9909.188034, 10657.05128, 10657.05128, 24866.45299, -24118.58974, -560.8974359, -27857.90598, 560.8974359, 560.8974359, -38327.99145, 26362.17949, -10657.05128, -36832.26496, 24118.58974, 24118.58974, -38327.99145, -10657.05128, 26362.17949 },
				{ -26362.17949, 81704.05983, 26362.17949, 560.8974359, -27857.90598, -560.8974359, -10657.05128, -38327.99145, -26362.17949, -24118.58974, 24866.45299, 560.8974359, -560.8974359, 24866.45299, 24118.58974, 26362.17949, -38327.99145, 10657.05128, 24118.58974, -36832.26496, -24118.58974, 10657.05128, 9909.188034, -10657.05128 },
				{ -26362.17949, 26362.17949, 81704.05983, -24118.58974, 560.8974359, 24866.45299, -10657.05128, -26362.17949, -38327.99145, 560.8974359, -560.8974359, -27857.90598, -560.8974359, 24118.58974, 24866.45299, 10657.05128, -10657.05128, 9909.188034, 24118.58974, -24118.58974, -36832.26496, 26362.17949, 10657.05128, -38327.99145 },
				{ 24866.45299, 560.8974359, -24118.58974, 81704.05983, 26362.17949, -26362.17949, 24866.45299, 24118.58974, -560.8974359, 9909.188034, -10657.05128, 10657.05128, -38327.99145, -26362.17949, -10657.05128, -27857.90598, -560.8974359, 560.8974359, -38327.99145, 10657.05128, 26362.17949, -36832.26496, -24118.58974, 24118.58974 },
				{ -560.8974359, -27857.90598, 560.8974359, 26362.17949, 81704.05983, -26362.17949, 24118.58974, 24866.45299, -560.8974359, 10657.05128, -38327.99145, 26362.17949, -26362.17949, -38327.99145, -10657.05128, 560.8974359, 24866.45299, -24118.58974, -10657.05128, 9909.188034, 10657.05128, -24118.58974, -36832.26496, 24118.58974 },
				{ -24118.58974, -560.8974359, 24866.45299, -26362.17949, -26362.17949, 81704.05983, 560.8974359, 560.8974359, -27857.90598, -10657.05128, 26362.17949, -38327.99145, 10657.05128, 10657.05128, 9909.188034, -560.8974359, -24118.58974, 24866.45299, 26362.17949, -10657.05128, -38327.99145, 24118.58974, 24118.58974, -36832.26496 },
				{ 9909.188034, -10657.05128, -10657.05128, 24866.45299, 24118.58974, 560.8974359, 81704.05983, 26362.17949, 26362.17949, 24866.45299, 560.8974359, 24118.58974, -36832.26496, -24118.58974, -24118.58974, -38327.99145, 10657.05128, -26362.17949, -27857.90598, -560.8974359, -560.8974359, -38327.99145, -26362.17949, 10657.05128 },
				{ 10657.05128, -38327.99145, -26362.17949, 24118.58974, 24866.45299, 560.8974359, 26362.17949, 81704.05983, 26362.17949, -560.8974359, -27857.90598, -560.8974359, -24118.58974, -36832.26496, -24118.58974, -10657.05128, 9909.188034, -10657.05128, 560.8974359, 24866.45299, 24118.58974, -26362.17949, -38327.99145, 10657.05128 },
				{ 10657.05128, -26362.17949, -38327.99145, -560.8974359, -560.8974359, -27857.90598, 26362.17949, 26362.17949, 81704.05983, 24118.58974, 560.8974359, 24866.45299, -24118.58974, -24118.58974, -36832.26496, -26362.17949, 10657.05128, -38327.99145, 560.8974359, 24118.58974, 24866.45299, -10657.05128, -10657.05128, 9909.188034 },
				{ 24866.45299, -24118.58974, 560.8974359, 9909.188034, 10657.05128, -10657.05128, 24866.45299, -560.8974359, 24118.58974, 81704.05983, -26362.17949, 26362.17949, -38327.99145, -10657.05128, -26362.17949, -36832.26496, 24118.58974, -24118.58974, -38327.99145, 26362.17949, 10657.05128, -27857.90598, 560.8974359, -560.8974359 },
				{ -24118.58974, 24866.45299, -560.8974359, -10657.05128, -38327.99145, 26362.17949, 560.8974359, -27857.90598, 560.8974359, -26362.17949, 81704.05983, -26362.17949, 10657.05128, 9909.188034, 10657.05128, 24118.58974, -36832.26496, 24118.58974, 26362.17949, -38327.99145, -10657.05128, -560.8974359, 24866.45299, -24118.58974 },
				{ -560.8974359, 560.8974359, -27857.90598, 10657.05128, 26362.17949, -38327.99145, 24118.58974, -560.8974359, 24866.45299, 26362.17949, -26362.17949, 81704.05983, -26362.17949, -10657.05128, -38327.99145, -24118.58974, 24118.58974, -36832.26496, -10657.05128, 10657.05128, 9909.188034, 560.8974359, -24118.58974, 24866.45299 },
				{ -27857.90598, -560.8974359, -560.8974359, -38327.99145, -26362.17949, 10657.05128, -36832.26496, -24118.58974, -24118.58974, -38327.99145, 10657.05128, -26362.17949, 81704.05983, 26362.17949, 26362.17949, 24866.45299, 560.8974359, 24118.58974, 9909.188034, -10657.05128, -10657.05128, 24866.45299, 24118.58974, 560.8974359 },
				{ 560.8974359, 24866.45299, 24118.58974, -26362.17949, -38327.99145, 10657.05128, -24118.58974, -36832.26496, -24118.58974, -10657.05128, 9909.188034, -10657.05128, 26362.17949, 81704.05983, 26362.17949, -560.8974359, -27857.90598, -560.8974359, 10657.05128, -38327.99145, -26362.17949, 24118.58974, 24866.45299, 560.8974359 },
				{ 560.8974359, 24118.58974, 24866.45299, -10657.05128, -10657.05128, 9909.188034, -24118.58974, -24118.58974, -36832.26496, -26362.17949, 10657.05128, -38327.99145, 26362.17949, 26362.17949, 81704.05983, 24118.58974, 560.8974359, 24866.45299, 10657.05128, -26362.17949, -38327.99145, -560.8974359, -560.8974359, -27857.90598 },
				{ -38327.99145, 26362.17949, 10657.05128, -27857.90598, 560.8974359, -560.8974359, -38327.99145, -10657.05128, -26362.17949, -36832.26496, 24118.58974, -24118.58974, 24866.45299, -560.8974359, 24118.58974, 81704.05983, -26362.17949, 26362.17949, 24866.45299, -24118.58974, 560.8974359, 9909.188034, 10657.05128, -10657.05128 },
				{ 26362.17949, -38327.99145, -10657.05128, -560.8974359, 24866.45299, -24118.58974, 10657.05128, 9909.188034, 10657.05128, 24118.58974, -36832.26496, 24118.58974, 560.8974359, -27857.90598, 560.8974359, -26362.17949, 81704.05983, -26362.17949, -24118.58974, 24866.45299, -560.8974359, -10657.05128, -38327.99145, 26362.17949 },
				{ -10657.05128, 10657.05128, 9909.188034, 560.8974359, -24118.58974, 24866.45299, -26362.17949, -10657.05128, -38327.99145, -24118.58974, 24118.58974, -36832.26496, 24118.58974, -560.8974359, 24866.45299, 26362.17949, -26362.17949, 81704.05983, -560.8974359, 560.8974359, -27857.90598, 10657.05128, 26362.17949, -38327.99145 },
				{ -36832.26496, 24118.58974, 24118.58974, -38327.99145, -10657.05128, 26362.17949, -27857.90598, 560.8974359, 560.8974359, -38327.99145, 26362.17949, -10657.05128, 9909.188034, 10657.05128, 10657.05128, 24866.45299, -24118.58974, -560.8974359, 81704.05983, -26362.17949, -26362.17949, 24866.45299, -560.8974359, -24118.58974 },
				{ 24118.58974, -36832.26496, -24118.58974, 10657.05128, 9909.188034, -10657.05128, -560.8974359, 24866.45299, 24118.58974, 26362.17949, -38327.99145, 10657.05128, -10657.05128, -38327.99145, -26362.17949, -24118.58974, 24866.45299, 560.8974359, -26362.17949, 81704.05983, 26362.17949, 560.8974359, -27857.90598, -560.8974359 },
				{ 24118.58974, -24118.58974, -36832.26496, 26362.17949, 10657.05128, -38327.99145, -560.8974359, 24118.58974, 24866.45299, 10657.05128, -10657.05128, 9909.188034, -10657.05128, -26362.17949, -38327.99145, 560.8974359, -560.8974359, -27857.90598, -26362.17949, 26362.17949, 81704.05983, -24118.58974, 560.8974359, 24866.45299 },
				{ -38327.99145, 10657.05128, 26362.17949, -36832.26496, -24118.58974, 24118.58974, -38327.99145, -26362.17949, -10657.05128, -27857.90598, -560.8974359, 560.8974359, 24866.45299, 24118.58974, -560.8974359, 9909.188034, -10657.05128, 10657.05128, 24866.45299, 560.8974359, -24118.58974, 81704.05983, 26362.17949, -26362.17949 },
				{ -10657.05128, 9909.188034, 10657.05128, -24118.58974, -36832.26496, 24118.58974, -26362.17949, -38327.99145, -10657.05128, 560.8974359, 24866.45299, -24118.58974, 24118.58974, 24866.45299, -560.8974359, 10657.05128, -38327.99145, 26362.17949, -560.8974359, -27857.90598, 560.8974359, 26362.17949, 81704.05983, -26362.17949 },
				{ 26362.17949, -10657.05128, -38327.99145, 24118.58974, 24118.58974, -36832.26496, 10657.05128, 10657.05128, 9909.188034, -560.8974359, -24118.58974, 24866.45299, 560.8974359, 560.8974359, -27857.90598, -10657.05128, 26362.17949, -38327.99145, -24118.58974, -560.8974359, 24866.45299, -26362.17949, -26362.17949, 81704.05983 }
			}));

			//var element = new Element() { ID = 0, ElementType = new Hexa8(Material1) };
			////var element = new Element() { ID = 0, ElementType = new Hexa8Fixed(Material1) };
			//foreach (var node in AbaqusToMSolveHexa8Nodes(NodeSet1)) element.AddNode(node);
			//writer.WriteToFile(element.ElementType.StiffnessMatrix(element), output, true);

			Assert.True(computedK.Equals(expectedK, 1e-8));
		}
	}
}
