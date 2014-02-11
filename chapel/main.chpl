use BlockDist;

// Define
const BND_LEFT: int = 1, BND_RIGHT: int = 2, BND_BOTTOM: int = 3, BND_TOP: int = 4;

// Structure
record problem {
    /* Problem variables */
    var alpha: real;
    var length_x, length_y, simulation_time: real;
    var nb_cell_x, nb_cell_y, nb_timestep: int;
    var dx, dy, dt: real;
    // Initial temperature
    var temp_initial: real;
    // Boundaries
    var bnd_type: [1..4] int;
    var bnd_value: [1..4] real;
    // Domain
    var physicalDomain, completeDomain: domain(2);

    /* Procedure: Reads problem configuration file
            f: file name
            */
    proc read_config_file(f: string) {
        // file is opened in read mode
        var configFile = open(f, iomode.r);
        var reader = configFile.reader();
        // TODO: modify the way the file is read using read() instead of readline or readln
        var line: string;
        // lecture du fichier
        reader.readline(line);
        alpha = reader.read(real);
        reader.readline(line);
        reader.readline(line);
        length_x = reader.read(real);
        nb_cell_x = reader.read(int);
        length_y = reader.read(real);
        nb_cell_y = reader.read(int);
        reader.readline(line);
        reader.readline(line);
        simulation_time = reader.read(real);
        nb_timestep = reader.read(int);
        reader.readline(line);
        reader.readline(line);
        temp_initial = reader.read(real);
        reader.readline(line);
        reader.readline(line);
        bnd_type(BND_LEFT) = reader.read(int);
        bnd_value(BND_LEFT) = reader.read(real);
        reader.readline(line);
        reader.readline(line);
        bnd_type(BND_RIGHT) = reader.read(int); 
        bnd_value(BND_RIGHT) = reader.read(real);
        reader.readline(line);
        reader.readline(line);
        bnd_type(BND_BOTTOM) = reader.read(int); 
        bnd_value(BND_BOTTOM) = reader.read(real);
        reader.readline(line);
        reader.readline(line);
        bnd_type(BND_TOP) = reader.read(int);
        bnd_value(BND_TOP) = reader.read(real);
        
        // file is closed
        reader.close();
        configFile.close();
        
        // based on read parameters, some characteristics are calculated
        dx = length_x/nb_cell_x;
        dy = length_y/nb_cell_y;
        dt = simulation_time/nb_timestep;
        // domain creation
        var domainTemp: domain(2) dmapped Block({1..nb_cell_x, 1..nb_cell_y}) = {1..nb_cell_x, 1..nb_cell_y};
        physicalDomain = domainTemp;
        completeDomain = physicalDomain.expand(1);
    }
    
    
    /* Procedure: Writes the case file for Paraview
            */
    proc write_case() {
        // file creation
        var caseFile = open("solut/ensight.case", iomode.cw);
        var w = caseFile.writer();

        w.writeln("FORMAT\ntype: ensight gold\n");
        w.writeln("GEOMETRY\nmodel: ensight.geo\n");
        w.writeln("VARIABLE\nscalar per element: Temperature temp_*******.ensight\n");
        w.writeln("TIME\ntime set: 1\nnumber of steps: ", nb_timestep, "\nfilename start number: 0\nfilename increment: 1\ntime values:");

        for i in 1..nb_timestep {
            w.writeln(" ", i * dt, new iostyle(realfmt=2, precision=5));
        }
        
        w.close();
        caseFile.close();
    }

    /* Procedure: Writes the geometry file for Paraview
                   */
    proc write_geo() {
        var geoFile = open("solut/ensight.geo", iomode.cw);
        var w = geoFile.writer();

        w.writeln("Ensight Gold Format file");
        w.writeln("Heat equation example program");
        w.writeln("node id off\nelement id off");
        w.writeln("part\n1");
        w.writeln("Mesh\nblock");
        w.writeln(nb_cell_x + 1, " ", nb_cell_y + 1, " ", 1);            

        for j in 0..nb_cell_y {
            for i in 0..nb_cell_x {
                w.writeln(dx * i, new iostyle(realfmt=2, precision=5));
            }
        }

        for j in 0..nb_cell_y {
            for i in 0..nb_cell_x {
                w.writeln(dy * j, new iostyle(realfmt=2, precision=5));
            }
        }

        for j in 0..nb_cell_y {
            for i in 0..nb_cell_x {
                w.writeln(i*0, new iostyle(realfmt=2, precision=5));
            }
        }

        w.close();
        geoFile.close();
    }
    
    /* Procedure: Writes the data file of each iteration for Paraview
            temp: reference to the domain
            it: iteration number
            */
    proc write_data(ref temp: [completeDomain] real, it: int) {
        // data file name containing iteration number
        var fileName: string;
        fileName.write("solut/temp_");
        // TODO: find a better way
        fileName.write(it/1000000);
        fileName.write(it/100000);
        fileName.write(it/10000);
        fileName.write(it/1000);
        fileName.write(it/100);
        fileName.write(it/10);
        fileName.write(it%10);
        fileName.write(".ensight");
        
        var dataFile = open(fileName, iomode.cw);
        var w = dataFile.writer();

        w.writeln("Temperature at cells");
        w.writeln("part\n1");
        w.writeln("block");

        for j in 1..nb_cell_y{
            for i in 1..nb_cell_x {
                w.writeln(temp(i, j), new iostyle(realfmt=2, precision=5));
            }
        }

        w.close();
        dataFile.close();
    }
}

/**************************************************/
// config variable: configuration file name
config var f: string = "config.dat";
config var verbose: bool = true;
config var write_data: bool = true;
var pb: problem;

pb.read_config_file(f);
if(write_data) then pb.write_case();
if(write_data) then pb.write_geo();

// Domain initialization
var temp, temp_old: [pb.completeDomain] real = pb.temp_initial; 

/* CFL CONDITION */
if(pb.alpha*pb.dt/min(pb.dx,pb.dy) >= 0.25){
    if(verbose) then writeln("CFL is too high:");
} else{
    if(verbose) then writeln("CFL is OK");
}

/* Stencil definition */
const stencil = ((0,0), (1,0), (0,1), (-1,0), (0,-1));

const A = (1-2*pb.alpha*pb.dt*(pb.dy/pb.dx+pb.dx/pb.dy)/pb.dx/pb.dy, pb.alpha*pb.dt/pb.dx/pb.dx, pb.alpha*pb.dt/pb.dy/pb.dy, pb.alpha*pb.dt/pb.dx/pb.dx, pb.alpha*pb.dt/pb.dy/pb.dy);

/* Iteration */
/*   Cette boucle n'est pas parallèlisée car on ne veut pas créer des processus parallèls s'occupant chacun d'un timestep vu que l'on veut qu'ils se suivent les uns les autres */
for k in 1..pb.nb_timestep {
    if(verbose) then writeln("Timestep#: ", k);
    
    // Ghost cell update
    forall j in 1..pb.nb_cell_y {
        temp_old(0, j) = pb.bnd_type(BND_LEFT) * (2 * pb.bnd_value(BND_LEFT) - temp_old(1, j)) + (1 - pb.bnd_type(BND_LEFT)) * (temp_old(1, j) - pb.dy * pb.bnd_value(BND_LEFT));
        temp_old(pb.nb_cell_x+1, j) = pb.bnd_type(BND_RIGHT) * (2 * pb.bnd_value(BND_RIGHT) - temp_old(pb.nb_cell_x, j)) + (1 - pb.bnd_type(BND_RIGHT)) * (temp_old(pb.nb_cell_x, j) - pb.dy * pb.bnd_value(BND_RIGHT));
    }
    forall i in 1..pb.nb_cell_x {
        temp_old(i, 0) = pb.bnd_type(BND_BOTTOM) * (2 * pb.bnd_value(BND_BOTTOM) - temp_old(i, 1)) + (1 - pb.bnd_type(BND_BOTTOM)) * (temp_old(i, 1) - pb.dx * pb.bnd_value(BND_BOTTOM));
        temp_old(i, pb.nb_cell_y+1) = pb.bnd_type(BND_TOP) * (2 * pb.bnd_value(BND_TOP) - temp_old(i, pb.nb_cell_y)) + (1 - pb.bnd_type(BND_TOP)) * (temp_old(i, pb.nb_cell_y) - pb.dx * pb.bnd_value(BND_TOP));
    }

    // Parallel calculation of the temperature of the domain
    forall cell in pb.physicalDomain {
        temp(cell) = A(1)*temp_old(cell + stencil(1)) + A(2)*temp_old(cell + stencil(2)) + A(3)*temp_old(cell + stencil(3)) + A(4)*temp_old(cell + stencil(4)) + A(5)*temp_old(cell + stencil(5));
        // la reduction semble TRES lente
        //temp(cell) = (+ reduce [i in 1..5] A(i)*temp_old(cell + stencil(i)));
    }
    
    // Data storage for next iteration
    temp_old = temp;
    
    // Data file writing
    if(write_data) then pb.write_data(temp, k-1);
}

