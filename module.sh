
name=$1
namecap=$(echo "$name" |  tr '[:lower:]' '[:upper:]' )

modulecpp="$name.cpp"
modulehpp="$name.hpp"

#test the existance of the module
if [ -d "src/$name" ];
then
    echo "Module $name exists."
else
    #create pertinent files and directories
    mkdir src/$name
    touch src/$name/$modulecpp
    touch src/$name/$modulehpp

    #write header guards and include statements
    printf "#ifndef $namecap"_HPP"\n" >> src/$name/$modulehpp
    printf "#define $namecap"_HPP"\n\n\n\n" >> src/$name/$modulehpp
    printf "#endif" >> src/$name/$modulehpp

    printf "#include \"$modulehpp\" " >> src/$name/$modulecpp

    echo "Module $name has been created."
fi
