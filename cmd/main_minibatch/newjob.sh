#!/bin/bash
#SBATCH -J Lattigo
#SBATCH -p defq
#SBATCH -A chenjingwei
#SBATCH -N 1
#SBATCH -w node01
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=80G

go run main.go
