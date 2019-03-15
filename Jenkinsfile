pipeline {
  agent { docker { image 'ubuntu-doxygen' } }
  stages {
    stage('Generate documentation') {
      steps {
        sh 'rm -rf gh-pages'
        sh 'git clone -b gh-pages https://github.com/BlackHolePerturbationToolkit/GremlinEq.git gh-pages'
        sh '( cat Doxyfile; echo "OUTPUT_DIRECTORY=gh-pages"; echo "HTML_OUTPUT=doc"; echo "GENERATE_LATEX=NO"; echo "SHORT_NAMES=YES";) | doxygen -'
        dir(path: 'gh-pages') {
          sh 'git add .'
          sh 'git -c user.name="BHPT Jenkins" -c user.email="" commit -am "Update documentation" || true'
          sh 'git push'
        }
      }
    }
  }
}