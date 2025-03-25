//valve3-2-01.cpp : 定义控制台应用程序的入口点。
//修改input文件需要修改对应的AddVAlve代码!!!"/100.0"

#include "stdafx.h"
#include "epanet2_2.h"
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "cmath"
#include "time.h"

#define ValveNum 3
#define Minmum 16
#define POPSIZE 100
#define MAXGENS 2000
#define FITNESS 10000
#define PatternNum 17
#define CoefDemand 1.0

int debug = 0;
const int VARIABLE = 2 * ValveNum + 3;
const double g = 9.81;
int generation, i, j, Nnodes, Nlinks, ReservoirIndex[3], resultoutput = 0;//Nnodes=4246, Nlinks=4845;
double Quality, Demand, Pressure, Elevation, MaxPressure, MinPressure, BaseDemand, patternindex, reservoirhead[3], reservoirdemand[3];
double LinkStatus, Length, startdiameter[4846];//从startroughness[1]开始存储，startdiameter[1]
int *presultoutput;
FILE *output;
FILE *PressureDistribution;
FILE *HeightDistribution;
FILE *PressureAccumulation;
FILE *NodeInformation;
FILE *PatternValue;
FILE *ReservoirHead1;
FILE *ReservoirHead2;
FILE *ReservoirHead3;
FILE *LeakageCurve;
FILE *NodeCostCurve;
FILE *PenaltyCurve;
FILE *TotalCostCurve;
FILE *wtpcostCurve;
EN_Project region;
struct NodeHeight {
	int Index;
	int SecondaryPattern;
	int Floor;
	double FloorHeight;
	double RegionLoss;
	double LossFactor;
	double Height;
};
struct ValveLink {
	char LinkID[20];
	int Number;
};
struct genotype {
	int gene[VARIABLE];
	double fitness;
	int upper[VARIABLE];
	int lower[VARIABLE];
	double rfitness;
	double cfitness;
};
struct NodeHeight Information[5000];//从Information[1]开始读入数据
struct ValveLink ValveName[ValveNum];
struct genotype population[POPSIZE + 1];
struct genotype newpopulation[POPSIZE + 1];

void initialize(void);
void inherit(void);
int randval(int, int);
void evaluate(void);
double TotalCost(int *, double *, int);
void Reservoir(void);
void AddValve(int, char *, double);
double Max(void);
double Min(void);
void RecoverValve(void);
void keep_the_best(void);
void report(void);
void select(void);
void crossover(double);
void exchange(int, int);
void mutate(double);
void elitist(void);

int main()
{
	srand((unsigned)time(NULL));
	//打开并读取inp文件
	EN_createproject(&region);
	EN_open(region, "Huzhou02.inp", "Huzhou02.rpt", "");
	EN_getnodeindex(region, "Reservior1", &ReservoirIndex[0]);
	EN_getnodeindex(region, "Reservior2", &ReservoirIndex[1]);
	EN_getnodeindex(region, "Reservior3", &ReservoirIndex[2]);
	EN_getcount(region, EN_LINKCOUNT, &Nlinks);
	EN_getcount(region, EN_NODECOUNT, &Nnodes);
	for (i = 1; i <= 4845; i++)
		EN_getlinkvalue(region, i, EN_DIAMETER, &startdiameter[i]);
	//打开文件与读入数据
	if ((output = fopen("output.txt", "w")) == NULL)
		exit(1);
	if ((PressureAccumulation = fopen("PressureAccumulation.txt", "w")) == NULL)
		exit(1);
	if ((PressureDistribution = fopen("PressureDistribution.txt", "w")) == NULL)
		exit(1);
	if ((HeightDistribution = fopen("HeightDistribution.txt", "w")) == NULL)
		exit(1);
	if ((LeakageCurve = fopen("LeakageCurve.csv", "w")) == NULL)
		exit(1);
	if ((NodeCostCurve = fopen("NodeCostCurve.csv", "w")) == NULL)
		exit(1);
	if ((PenaltyCurve = fopen("PenaltyCurve.csv", "w")) == NULL)
		exit(1);
	if ((TotalCostCurve = fopen("TotalCostCurve.csv", "w")) == NULL)
		exit(1);
	if ((wtpcostCurve = fopen("wtpcostCurve.csv", "w")) == NULL)
		exit(1);
	if ((NodeInformation = fopen("NodeInformation.txt", "r")) == NULL) {
		fprintf(output, "\n Cannot open NodeInformation file! \n");
		exit(1);
	}
	for (i = 1; i <= Nnodes; i++) {//这里只循环4246次了，NI文件里还是4266个节点
		fscanf(NodeInformation, "%d", &Information[i].Index);
		fscanf(NodeInformation, "%d", &Information[i].SecondaryPattern);
		fscanf(NodeInformation, "%d", &Information[i].Floor);
		fscanf(NodeInformation, "%lf", &Information[i].FloorHeight);
		fscanf(NodeInformation, "%lf", &Information[i].RegionLoss);
		fscanf(NodeInformation, "%lf", &Information[i].LossFactor);
		fscanf(NodeInformation, "%lf", &Information[i].Height);
	}
	fclose(NodeInformation);
	for (i = 1; i <= Nnodes - 3; i++) {
		EN_getnodevalue(region, i, EN_INITQUAL, &Quality);
		if (Information[i].Height != Quality) {
			fprintf(output, "\n Error! Information.Height != Quality \n");
			exit(1);
		}
	}
	//种群初始化、评估、择优、报告
	presultoutput = &resultoutput;
	generation = 0;
	initialize();
	// inherit();
	fprintf(output, "generation, best_fitness, average, standard_deviation\n");
	printf("generation, best_fitness, average, standard_deviation\n");
	evaluate();
	keep_the_best();
	report();
	//进化
	double Pmutation, Pcrossover;
	for (generation = 1; generation <= MAXGENS; generation++) {
		if (generation < MAXGENS / 2 + 1) {
			Pmutation = 0.3 + 0.3*generation * 2 / MAXGENS;
			Pcrossover = 0.6 - 0.3*generation * 2 / MAXGENS;
		}
		else {
			Pmutation = 0.6;
			Pcrossover = 0.3;
		}
		select();
		crossover(Pcrossover);
		mutate(Pmutation);
		evaluate();
		elitist();
		report();
		for (i = 0; i < 3; i++) {
			fprintf(output, "Reservoirhead[%d] = %d\n", i, population[POPSIZE].gene[i]);
			printf("Reservoirhead[%d] = %d\n", i, population[POPSIZE].gene[i]);
		}
		for (i = 0; i < ValveNum; i++) {
			fprintf(output, "ValveID[%d] = %d\n", i, population[POPSIZE].gene[i + 3]);
			printf("ValveID[%d] = %d\n", i, population[POPSIZE].gene[i + 3]);
		}
		for (i = 0; i < ValveNum; i++) {
			fprintf(output, "vdiameter[%d] = %d\n", i, population[POPSIZE].gene[i + 3 + ValveNum]);
			printf("vdiameter[%d] = %d\n", i, population[POPSIZE].gene[i + 3 + ValveNum]);
		}
	}
	//打印模拟完成以及最优基因
	double error = 0;
	fprintf(output, "\n Simulation Completed! \n");
	printf("\n Simulation Completed! \n");
	*presultoutput = 1;
	int ValveID[ValveNum];
	double vdiameter[ValveNum];
	for (i = 0; i < 3; i++)
		reservoirhead[i] = population[POPSIZE].gene[i] / 10.0;
	for (i = 0; i < ValveNum; i++)
		ValveID[i] = population[POPSIZE].gene[i + 3];
	for (i = 0; i < ValveNum; i++)
		vdiameter[i] = population[POPSIZE].gene[i + 3 + ValveNum];
	error = FITNESS - TotalCost(ValveID, vdiameter, POPSIZE);
	if (error == 1) {
		printf("\n population[POPSIZE] error! \n");
		fprintf(output, "\n population[POPSIZE] error! \n");
	}
	int linkindex;
	for (i = 0; i < ValveNum; i++) {
		sprintf(ValveName[i].LinkID, "300MM%d", population[POPSIZE].gene[i + 3]);
		AddValve(i, ValveName[i].LinkID, vdiameter[i]);
		EN_getlinkindex(region, ValveName[i].LinkID, &linkindex);
		EN_setlinkvalue(region, linkindex, EN_KBULK, 2.0);
	}
	EN_saveinpfile(region, "Huzhou03.inp");
	//压力分布统计
	EN_solveH(region);
	int Accumulation[80];
	double AccumulationSum = 0, AccumulationPlot[80];
	for (j = 0; j < 80; j++)
		Accumulation[j] = 0;
	for (i = 1; i <= Nnodes - 3; i++) {
		EN_getnodevalue(region, i, EN_PRESSURE, &Pressure);
		for (j = 0; j < 80; j++)
			if (Pressure >= (15 + 0.5*j) && Pressure < (15.5 + 0.5*j))
				Accumulation[j] = Accumulation[j] + 1;
	}
	for (j = 0; j < 80; j++)
		AccumulationSum = AccumulationSum + Accumulation[j];
	AccumulationPlot[0] = Accumulation[0] / AccumulationSum;
	for (j = 1; j < 80; j++)
		AccumulationPlot[j] = AccumulationPlot[j - 1] + Accumulation[j] / AccumulationSum;
	for (j = 0; j < 80; j++)
		fprintf(PressureAccumulation, "%lf\n", AccumulationPlot[j]);
	for (i = 1; i <= Nnodes - 3; i++) {
		EN_getnodevalue(region, i, EN_PRESSURE, &Pressure);
		fprintf(PressureDistribution, "%lf\n", Pressure);
		EN_getnodevalue(region, i, EN_ELEVATION, &Elevation);
		fprintf(HeightDistribution, "%lf\n", Pressure + Elevation);
	}
	//关闭文件与项目
	fclose(HeightDistribution);
	fclose(PressureDistribution);
	fclose(PressureAccumulation);
	fclose(output);
	fclose(LeakageCurve);
	fclose(NodeCostCurve);
	fclose(PenaltyCurve);
	fclose(TotalCostCurve);
	fclose(wtpcostCurve);
	EN_closeH(region);
	EN_close(region);
	EN_deleteproject(region);

	return 0;
}
//依据input.txt初始化POPSIZE个种群中的VARIABLE个变量
void initialize(void) {
	for (i = 1; i <= Nnodes - 3; i++) {
		EN_getnodevalue(region, i, EN_BASEDEMAND, &BaseDemand);
		EN_setnodevalue(region, i, EN_BASEDEMAND, CoefDemand * BaseDemand);
	}
	FILE *input;
	int i, j, lbound, ubound;
	if ((input = fopen("input.txt", "r")) == NULL) {
		fprintf(output, "\n Cannot open input file! \n");
		exit(1);
	}
	for (i = 0; i < VARIABLE; i++) {
		fscanf(input, "%d", &lbound);
		fscanf(input, "%d", &ubound);
		for (j = 0; j < POPSIZE; j++) {
			population[j].fitness = 0;
			population[j].rfitness = 0;
			population[j].cfitness = 0;
			population[j].lower[i] = lbound;
			population[j].upper[i] = ubound;
			population[j].gene[i] = randval(population[j].lower[i], population[j].upper[i]);
		}
	}
	population[POPSIZE].fitness = 0;
	fclose(input);
}
//继承阀门数量-1的阀门配置
void inherit(void) {
	FILE *inherit;
	int last;
	if ((inherit = fopen("inherit.txt", "r")) == NULL) {
		fprintf(output, "\n Cannot open input file! \n");
		exit(1);
	}
	for (i = 0; i < 3; i++) {
		fscanf(inherit, "%d", &last);
		population[0].gene[i] = last;
	}
	for (i = 0; i < ValveNum - 1; i++) {
		fscanf(inherit, "%d", &last);
		population[0].gene[i + 3] = last;
		fscanf(inherit, "%d", &last);
		population[0].gene[i + 3 + ValveNum] = last;
	}
	population[0].gene[VARIABLE - 1] = 100;
	fclose(inherit);
}
//在两个值之间产生随机数并返回
int randval(int low, int high) {
	int val;
	val = rand() % (high - low + 1) + low;
	return(val);
}
//评估个体适应度
void evaluate(void) {
	int mem, ValveID[ValveNum];
	double vdiameter[ValveNum];
	for (mem = 0; mem < POPSIZE; mem++) {
		for (i = 0; i < 3; i++)
			reservoirhead[i] = population[mem].gene[i] / 10.0;
		for (i = 0; i < ValveNum; i++)
			ValveID[i] = population[mem].gene[i + 3];
		for (i = 0; i < ValveNum; i++)
			vdiameter[i] = population[mem].gene[i + 3 + ValveNum];
		population[mem].fitness = FITNESS - TotalCost(ValveID, vdiameter, mem);
		if (population[mem].fitness <= 0)
			population[mem].fitness = 1;
	}
}
//年度总费用
double TotalCost(int ValveID[ValveNum], double vdiameter[ValveNum], int k) {
	//设置水库水头
	Reservoir();
	//添加阀门
	for (i = 0; i < ValveNum; i++) {
		sprintf(ValveName[i].LinkID, "300MM%d", ValveID[i]);
		AddValve(i, ValveName[i].LinkID, vdiameter[i]);
	}
	EN_solveH(region);
	MaxPressure = Max();
	MinPressure = Min();
	if (MinPressure < 0) {
		printf("\n There is a negative pressure! \n");
		RecoverValve();
		return FITNESS - 1;
	}
	//延时模拟及漏损、二供费用、惩罚值计算，在该版本中去除延时模拟，进行单步分析
	int NodeA, NodeB;
	double Leakage = 0, PressureA, PressureB, NodeFloor, NodeCost = 0, Penalty = 0;
	//漏损费用计算
	for (i = 1; i <= Nlinks; i++) {
		EN_getlinkvalue(region, i, EN_INITSTATUS, &LinkStatus);
		if (LinkStatus == 1) {//如果管道是开启状态
			EN_getlinkvalue(region, i, EN_LENGTH, &Length);
			EN_getlinknodes(region, i, &NodeA, &NodeB);
			EN_getnodevalue(region, NodeA, EN_PRESSURE, &PressureA);
			EN_getnodevalue(region, NodeB, EN_PRESSURE, &PressureB);
			Leakage = Leakage + (pow((PressureA + PressureB) / 2, 1.18))*Length * 1.3 / 100000 * 3 * 24 * 365 / 10000;//Leakage:m3/h, P&L:m, 1.3/100000表示漏损系数.
		}
	}
	//二供费用、惩罚值计算
	int ppressure = 0, phouse = 0;
	for (i = 1; i <= Nnodes - 3; i++) {
		EN_getnodevalue(region, i, EN_INITQUAL, &Quality);
		EN_getnodevalue(region, i, EN_PRESSURE, &Pressure);
		EN_getnodevalue(region, i, EN_DEMAND, &Demand);
		if (Pressure < Quality && Demand > 0) {
			NodeFloor = floor((Pressure - 10 - Information[i].RegionLoss + Information[i].FloorHeight)
				/ (Information[i].FloorHeight + Information[i].LossFactor*0.3));//floor向下取整，依据节点需求水压计算公式计算水压对应层数，10表示1+5+4，对应公式z+s+hflowmeter
			if (Information[i].SecondaryPattern == 0) {//无水池有水箱
				if (Information[i].Floor <= 12)
					NodeCost = NodeCost + Demand*(Information[i].Floor - NodeFloor) / Information[i].Floor*(Quality - Pressure) / 0.9 / 0.9;//两个0.9分别表示泵组效率和电机效率
				if (Information[i].Floor > 12 && Information[i].Floor <= 20) {
					if (NodeFloor <= 12)
						NodeCost = NodeCost + Demand*(Information[i].Floor - 12) / Information[i].Floor*(Quality - Pressure) / 0.9 / 0.9
						+ Demand*(12 - NodeFloor) / Information[i].Floor*(11 * Information[i].FloorHeight + 10 + Information[i].RegionLoss + Information[i].LossFactor*0.3 * 12 - Pressure) / 0.9 / 0.9;
					if (NodeFloor > 12)
						NodeCost = NodeCost + Demand*(Information[i].Floor - NodeFloor) / Information[i].Floor*(Quality - Pressure) / 0.9 / 0.9;
				}
				if (Information[i].Floor > 20) {
					if (NodeFloor <= 12)
						NodeCost = NodeCost + Demand*(Information[i].Floor - 20) / Information[i].Floor*(Quality - Pressure) / 0.9 / 0.9
						+ Demand*(20 - 12) / Information[i].Floor*(19 * Information[i].FloorHeight + 10 + Information[i].RegionLoss + Information[i].LossFactor*0.3 * 20 - Pressure) / 0.9 / 0.9
						+ Demand*(12 - NodeFloor) / Information[i].Floor*(11 * Information[i].FloorHeight + 10 + Information[i].RegionLoss + Information[i].LossFactor*0.3 * 12 - Pressure) / 0.9 / 0.9;
					if (NodeFloor > 12 && NodeFloor <= 20)
						NodeCost = NodeCost + Demand*(Information[i].Floor - 20) / Information[i].Floor*(Quality - Pressure) / 0.9 / 0.9
						+ Demand*(20 - NodeFloor) / Information[i].Floor*(19 * Information[i].FloorHeight + 10 + Information[i].RegionLoss + Information[i].LossFactor*0.3 * 20 - Pressure) / 0.9 / 0.9;
					if (NodeFloor > 20)
						NodeCost = NodeCost + Demand*(Information[i].Floor - NodeFloor) / Information[i].Floor*(Quality - Pressure) / 0.9 / 0.9;
				}
			}
			if (Information[i].SecondaryPattern == 1) {//无水池无水箱
				if (Information[i].Floor <= 12)
					NodeCost = NodeCost + Demand*(Information[i].Floor - NodeFloor) / Information[i].Floor*(Quality - Pressure) / 0.9 / 0.45;
				if (Information[i].Floor > 12 && Information[i].Floor <= 20) {
					if (NodeFloor <= 12)
						NodeCost = NodeCost + Demand*(Information[i].Floor - 12) / Information[i].Floor*(Quality - Pressure) / 0.9 / 0.45
						+ Demand*(12 - NodeFloor) / Information[i].Floor*(11 * Information[i].FloorHeight + 10 + Information[i].RegionLoss + Information[i].LossFactor*0.3 * 12 - Pressure) / 0.9 / 0.45;
					if (NodeFloor > 12)
						NodeCost = NodeCost + Demand*(Information[i].Floor - NodeFloor) / Information[i].Floor*(Quality - Pressure) / 0.9 / 0.45;
				}
				if (Information[i].Floor > 20) {
					if (NodeFloor <= 12)
						NodeCost = NodeCost + Demand*(Information[i].Floor - 20) / Information[i].Floor*(Quality - Pressure) / 0.9 / 0.45
						+ Demand*(20 - 12) / Information[i].Floor*(19 * Information[i].FloorHeight + 10 + Information[i].RegionLoss + Information[i].LossFactor*0.3 * 20 - Pressure) / 0.9 / 0.45
						+ Demand*(12 - NodeFloor) / Information[i].Floor*(11 * Information[i].FloorHeight + 10 + Information[i].RegionLoss + Information[i].LossFactor*0.3 * 12 - Pressure) / 0.9 / 0.45;
					if (NodeFloor > 12 && NodeFloor <= 20)
						NodeCost = NodeCost + Demand*(Information[i].Floor - 20) / Information[i].Floor*(Quality - Pressure) / 0.9 / 0.45
						+ Demand*(20 - NodeFloor) / Information[i].Floor*(19 * Information[i].FloorHeight + 10 + Information[i].RegionLoss + Information[i].LossFactor*0.3 * 20 - Pressure) / 0.9 / 0.45;
					if (NodeFloor > 20)
						NodeCost = NodeCost + Demand*(Information[i].Floor - NodeFloor) / Information[i].Floor*(Quality - Pressure) / 0.9 / 0.45;
				}
			}
			if (Information[i].SecondaryPattern == 2) {//有水池有水箱
				if (Information[i].Floor <= 12)
					NodeCost = NodeCost + Demand*(Information[i].Floor - NodeFloor) / Information[i].Floor*Quality / 0.9 / 0.9;
				if (Information[i].Floor > 12 && Information[i].Floor <= 20) {
					if (NodeFloor <= 12)
						NodeCost = NodeCost + Demand*(Information[i].Floor - 12) / Information[i].Floor*Quality / 0.9 / 0.9
						+ Demand*(12 - NodeFloor) / Information[i].Floor*(11 * Information[i].FloorHeight + 10 + Information[i].RegionLoss + Information[i].LossFactor*0.3 * 12) / 0.9 / 0.9;
					if (NodeFloor > 12)
						NodeCost = NodeCost + Demand*(Information[i].Floor - NodeFloor) / Information[i].Floor*Quality / 0.9 / 0.9;
				}
				if (Information[i].Floor > 20) {
					if (NodeFloor <= 12)
						NodeCost = NodeCost + Demand*(Information[i].Floor - 20) / Information[i].Floor*Quality / 0.9 / 0.9
						+ Demand*(20 - 12) / Information[i].Floor*(19 * Information[i].FloorHeight + 10 + Information[i].RegionLoss + Information[i].LossFactor*0.3 * 20) / 0.9 / 0.9
						+ Demand*(12 - NodeFloor) / Information[i].Floor*(11 * Information[i].FloorHeight + 10 + Information[i].RegionLoss + Information[i].LossFactor*0.3 * 12) / 0.9 / 0.9;
					if (NodeFloor > 12 && NodeFloor <= 20)
						NodeCost = NodeCost + Demand*(Information[i].Floor - 20) / Information[i].Floor*Quality / 0.9 / 0.9
						+ Demand*(20 - NodeFloor) / Information[i].Floor*(19 * Information[i].FloorHeight + 10 + Information[i].RegionLoss + Information[i].LossFactor*0.3 * 20) / 0.9 / 0.9;
					if (NodeFloor > 20)
						NodeCost = NodeCost + Demand*(Information[i].Floor - NodeFloor) / Information[i].Floor*Quality / 0.9 / 0.9;
				}
			}
			if (Information[i].SecondaryPattern == 3) {//有水池无水箱
				if (Information[i].Floor <= 12)
					NodeCost = NodeCost + Demand*(Information[i].Floor - NodeFloor) / Information[i].Floor*Quality / 0.9 / 0.45;
				if (Information[i].Floor > 12 && Information[i].Floor <= 20) {
					if (NodeFloor <= 12)
						NodeCost = NodeCost + Demand*(Information[i].Floor - 12) / Information[i].Floor*Quality / 0.9 / 0.45
						+ Demand*(12 - NodeFloor) / Information[i].Floor*(11 * Information[i].FloorHeight + 10 + Information[i].RegionLoss + Information[i].LossFactor*0.3 * 12) / 0.9 / 0.45;
					if (NodeFloor > 12)
						NodeCost = NodeCost + Demand*(Information[i].Floor - NodeFloor) / Information[i].Floor*Quality / 0.9 / 0.45;
				}
				if (Information[i].Floor > 20) {
					if (NodeFloor <= 12)
						NodeCost = NodeCost + Demand*(Information[i].Floor - 20) / Information[i].Floor*Quality / 0.9 / 0.45
						+ Demand*(20 - 12) / Information[i].Floor*(19 * Information[i].FloorHeight + 10 + Information[i].RegionLoss + Information[i].LossFactor*0.3 * 20) / 0.9 / 0.45
						+ Demand*(12 - NodeFloor) / Information[i].Floor*(11 * Information[i].FloorHeight + 10 + Information[i].RegionLoss + Information[i].LossFactor*0.3 * 12) / 0.9 / 0.45;
					if (NodeFloor > 12 && NodeFloor <= 20)
						NodeCost = NodeCost + Demand*(Information[i].Floor - 20) / Information[i].Floor*Quality / 0.9 / 0.45
						+ Demand*(20 - NodeFloor) / Information[i].Floor*(19 * Information[i].FloorHeight + 10 + Information[i].RegionLoss + Information[i].LossFactor*0.3 * 20) / 0.9 / 0.45;
					if (NodeFloor > 20)
						NodeCost = NodeCost + Demand*(Information[i].Floor - NodeFloor) / Information[i].Floor*Quality / 0.9 / 0.45;
				}
			}
		}
		if (Pressure < Minmum) {
			Penalty = Penalty + (Minmum - Pressure) * 24;
			ppressure = ppressure + 1;
		}
		if (ppressure > 5) {
			printf("ppressure error! \n");
			RecoverValve();
			return FITNESS - 1;
		}
		if (Pressure < Quality && Demand > 0) {//市政直供建筑
			if (Information[i].FloorHeight == 3.0 && Information[i].Floor <= 6) {
				Penalty = Penalty + (Quality - Pressure) * 24;
				phouse = phouse + 1;
			}
			if (Information[i].FloorHeight == 4.2 && Information[i].Floor <= 4) {
				Penalty = Penalty + (Quality - Pressure) * 24;
				phouse = phouse + 1;
			}
		}
		if (phouse > 25) {
			printf("phouse error! \n");
			RecoverValve();
			return FITNESS - 1;
		}
	}
	NodeCost = NodeCost / 1000 * 1000 * g / 1000 * 0.6 * 24 * 365 / 10000;
	//减压阀费用计算
	int valveindex;
	double ValveCost = 0;
	for (i = 0; i < ValveNum; i++) {
		EN_getlinkindex(region, ValveName[i].LinkID, &valveindex);
		if (startdiameter[valveindex] <= 350)
			ValveCost = ValveCost + 2.66;
		else if (startdiameter[valveindex] <= 500)
			ValveCost = ValveCost + 4.32;
		else if (startdiameter[valveindex] <= 700)
			ValveCost = ValveCost + 6.175;
		else if (startdiameter[valveindex] <= 900)
			ValveCost = ValveCost + 9.34;
		else if (startdiameter[valveindex] <= 1100)
			ValveCost = ValveCost + 16.43;
		else if (startdiameter[valveindex] <= 1300)
			ValveCost = ValveCost + 26.52;
		else if (startdiameter[valveindex] <= 1600)
			ValveCost = ValveCost + 36.87;
		else
			ValveCost = ValveCost + 79.86;
	}
	ValveCost = ValveCost * 1.1 / 20;//预期使用年限25年，年维护费用系数0.1
									 //水厂费用
	double ReservoirStart[3], wtpcost = 0, ReserviorBeyond[3];
	ReservoirStart[0] = 571.64*CoefDemand;
	ReservoirStart[1] = 453.14*CoefDemand;
	ReservoirStart[2] = 1897.44*CoefDemand;
	for (i = 0; i < 3; i++)
		EN_getnodevalue(region, ReservoirIndex[i], EN_DEMAND, &reservoirdemand[i]);
	for (i = 0; i < 3; i++) {
		if (fabs(reservoirdemand[i]) < 0.9*ReservoirStart[i] && fabs(reservoirdemand[i]) >= 0.8*ReservoirStart[i])
		{
			ReserviorBeyond[i] = 664 * pow(((0.9*ReservoirStart[i] - fabs(reservoirdemand[i])) / 1000 * 8.64), 0.8) / 20
				+ 542.75 * pow(((0.9*ReservoirStart[i] - fabs(reservoirdemand[i])) / 1000), 0.41) / 20;
			ReserviorBeyond[i] = ReserviorBeyond[i] * 2 + 100;
		}
		if (fabs(reservoirdemand[i]) > 1.1*ReservoirStart[i] && fabs(reservoirdemand[i]) <= 1.2*ReservoirStart[i])
		{
			ReserviorBeyond[i] = 664 * pow(((fabs(reservoirdemand[i]) - 1.1*ReservoirStart[i]) / 1000 * 8.64), 0.8) / 20
				+ 542.75 * pow(((fabs(reservoirdemand[i]) - 1.1*ReservoirStart[i]) / 1000), 0.41) / 20;
			ReserviorBeyond[i] = ReserviorBeyond[i] * 2 + 100;
		}
		if (fabs(reservoirdemand[i]) < 0.8*ReservoirStart[i])
		{
			ReserviorBeyond[i] = 664 * pow(((0.9*ReservoirStart[i] - fabs(reservoirdemand[i])) / 1000 * 8.64), 0.8) / 20
				+ 542.75 * pow(((0.9*ReservoirStart[i] - fabs(reservoirdemand[i])) / 1000), 0.41) / 20;
			ReserviorBeyond[i] = ReserviorBeyond[i] + 100;
		}
		if (fabs(reservoirdemand[i]) > 1.2*ReservoirStart[i])
		{
			ReserviorBeyond[i] = 664 * pow(((fabs(reservoirdemand[i]) - 1.1*ReservoirStart[i]) / 1000 * 8.64), 0.8) / 20
				+ 542.75 * pow(((fabs(reservoirdemand[i]) - 1.1*ReservoirStart[i]) / 1000), 0.41) / 20;
			ReserviorBeyond[i] = ReserviorBeyond[i] + 100;
		}
		if (fabs(reservoirdemand[i]) >= 0.9*ReservoirStart[i] && fabs(reservoirdemand[i]) <= 1.1*ReservoirStart[i])
			ReserviorBeyond[i] = 0;
		wtpcost = wtpcost + fabs(reservoirdemand[i]) * reservoirhead[i] / 1000 / 0.72 * 1000 * g / 1000 * 0.6 * 24 * 365 / 10000 + ReserviorBeyond[i];
	}
	//总成本计算
	double CostTotal;
	CostTotal = Leakage + NodeCost + Penalty + ValveCost + wtpcost;
	//if (*presultoutput != 0)
	fprintf(output, "population[%d]: Leakage = %lf, NodeCost = %lf, Penalty = %lf, ValveCost = %lf, CostTotal = %lf, wtpcost = %lf\n", k, Leakage, NodeCost, Penalty, ValveCost, CostTotal, wtpcost);
	printf("population[%d]: Leakage = %lf, NodeCost = %lf, Penalty = %lf, ValveCost = %lf, CostTotal = %lf, wtpcost = %lf\n", k, Leakage, NodeCost, Penalty, ValveCost, CostTotal, wtpcost);
	if (*presultoutput != 0) {
		fprintf(LeakageCurve, "%d,%lf\n", generation, Leakage);
		fprintf(NodeCostCurve, "%d,%lf\n", generation, NodeCost);
		fprintf(PenaltyCurve, "%d,%lf\n", generation, Penalty);
		fprintf(TotalCostCurve, "%d,%lf\n", generation, CostTotal);
		fprintf(wtpcostCurve, "%d,%lf\n", generation, wtpcost);
	}
	RecoverValve();
	return CostTotal;
}
//设置水库水头
void Reservoir(void) {
	EN_setnodevalue(region, 4244, EN_ELEVATION, reservoirhead[0] + 4.0);
	EN_setnodevalue(region, 4245, EN_ELEVATION, reservoirhead[1] + 5.0);
	EN_setnodevalue(region, 4246, EN_ELEVATION, reservoirhead[2] + 3.0);
}
//添加阀门
void AddValve(int ValveNumber, char ValveID[20], double vdiameter) {//此处第二个参数ValveID不需要传递，其为全局变量
	int linkindex;
	double Diameter;
	EN_getlinkindex(region, ValveID, &linkindex);
	EN_getlinkvalue(region, linkindex, EN_DIAMETER, &Diameter);
	EN_setlinkvalue(region, linkindex, EN_DIAMETER, vdiameter / 100.0 * Diameter);
}
//返回最大节点压力
double Max(void) {
	MaxPressure = -10;
	for (i = 1; i <= Nnodes - 3; i++) {
		EN_getnodevalue(region, i, EN_PRESSURE, &Pressure);
		if (Pressure > MaxPressure)
			MaxPressure = Pressure;
	}
	return MaxPressure;
}
//返回最小节点压力
double Min(void) {
	MinPressure = 1000;
	for (i = 1; i <= Nnodes - 3; i++) {
		EN_getnodevalue(region, i, EN_PRESSURE, &Pressure);
		if (Pressure < MinPressure)
			MinPressure = Pressure;
	}
	return MinPressure;
}
//回收阀门
void RecoverValve(void) {
	int linkindex;
	for (i = 0; i < ValveNum; i++) {
		EN_getlinkindex(region, ValveName[i].LinkID, &linkindex);
		EN_setlinkvalue(region, linkindex, EN_DIAMETER, startdiameter[linkindex]);
	}
}
//存储最优基因及适应度
void keep_the_best(void) {
	int mem, best;
	for (mem = 0; mem < POPSIZE; mem++) {
		if (population[mem].fitness > population[POPSIZE].fitness) {
			best = mem;
			population[POPSIZE].fitness = population[mem].fitness;
		}
	}
	for (i = 0; i < VARIABLE; i++)
		population[POPSIZE].gene[i] = population[best].gene[i];
}
//报告种群统计数据
void report(void) {
	double sum = 0, sum_square = 0, average, square_sum, standard_deviation;
	for (i = 0; i < POPSIZE; i++) {
		sum = sum + population[i].fitness;
		sum_square = sum_square + population[i].fitness*population[i].fitness;//平方和
	}
	average = sum / POPSIZE;
	square_sum = average*average*POPSIZE;//和的平方除以POPSIZE
	standard_deviation = sqrt((sum_square - square_sum) / (POPSIZE - 1));//样本标准差
	fprintf(output, "%d, %f, %f, %f\n", generation, population[POPSIZE].fitness, average, standard_deviation);
	printf("%d, %f, %f, %f\n", generation, population[POPSIZE].fitness, average, standard_deviation);
	//统计每代最优解的各项费用
	*presultoutput = 1;
	int ValveID[ValveNum];
	double vdiameter[ValveNum], error;
	for (i = 0; i < 3; i++)
		reservoirhead[i] = population[POPSIZE].gene[i] / 10.0;
	for (i = 0; i < ValveNum; i++)
		ValveID[i] = population[POPSIZE].gene[i + 3];
	for (i = 0; i < ValveNum; i++)
		vdiameter[i] = population[POPSIZE].gene[i + 3 + ValveNum];
	error = FITNESS - TotalCost(ValveID, vdiameter, POPSIZE);
	if (error == 1) {
		printf("\n population[POPSIZE] error! \n");
		fprintf(output, "\n population[POPSIZE] error! \n");
	}
	*presultoutput = 0;
}
//轮盘赌选择个体作为父代
void select(void) {
	int mem;
	double sum = 0, p;
	for (mem = 0; mem < POPSIZE; mem++)
		sum = sum + population[mem].fitness;
	for (mem = 0; mem < POPSIZE; mem++)
		population[mem].rfitness = population[mem].fitness / sum;
	population[0].cfitness = population[0].rfitness;
	for (mem = 1; mem < POPSIZE; mem++)
		population[mem].cfitness = population[mem - 1].cfitness + population[mem].rfitness;
	for (i = 0; i < POPSIZE; i++) {
		p = rand() % 1001 / 1000.0;
		if (p < population[0].cfitness)
			newpopulation[i] = population[0];
		else
			for (j = 1; j < POPSIZE; j++)
				if (p >= population[j - 1].cfitness && p < population[j].cfitness)
					newpopulation[i] = population[j];
	}
	for (i = 0; i < POPSIZE; i++)
		population[i] = newpopulation[i];
}
//交叉部分基因
void crossover(double probability) {
	int mem, flag = 1, first;
	double x;
	for (mem = 0; mem < POPSIZE; mem++) {
		x = rand() % 1001 / 1000.0;
		if (x < probability) {
			if (flag % 2 == 0)
				exchange(first, mem);
			else
				first = mem;
			flag = flag + 1;
		}
	}
}
//实现交换功能，仅交换前v-1个变量，随机确定个数
void exchange(int first, int second) {
	int point, temp;
	point = (rand() % (VARIABLE - 1)) + 1;
	for (i = 0; i < point; i++) {
		temp = population[first].gene[i];
		population[first].gene[i] = population[second].gene[i];
		population[second].gene[i] = temp;
	}
}
//一定概率重新随机基因
void mutate(double probability) {
	double x;
	for (i = 0; i < POPSIZE; i++) {//遍历每一个个体的每一个基因
		for (j = 0; j < VARIABLE; j++) {
			x = rand() % 1001 / 1000.0;
			if (x < probability)
				population[i].gene[j] = randval(population[i].lower[j], population[i].upper[j]);
		}
	}
}
//精英化
void elitist(void) {
	int bestmem = 0, worstmem = 0;
	double best = population[0].fitness, worst = population[0].fitness;
	for (i = 0; i < POPSIZE; i++) {
		if (population[i].fitness > best) {
			best = population[i].fitness;
			bestmem = i;
		}
		if (population[i].fitness < worst) {
			worst = population[i].fitness;
			worstmem = i;
		}
	}
	if (best > population[POPSIZE].fitness) {
		for (i = 0; i < VARIABLE; i++)
			population[POPSIZE].gene[i] = population[bestmem].gene[i];
		population[POPSIZE].fitness = population[bestmem].fitness;
	}
	else {
		for (i = 0; i < VARIABLE; i++)
			population[worstmem].gene[i] = population[POPSIZE].gene[i];
		population[worstmem].fitness = population[POPSIZE].fitness;
	}
}