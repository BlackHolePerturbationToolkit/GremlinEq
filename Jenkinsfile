pipeline {
  agent { docker { image 'ubuntu-doxygen' } }
  stages {
    stage('Generate documentation') {
      steps {
        withCredentials([sshUserPrivateKey(credentialsId: 'gremlin', keyFileVariable: 'SSH_KEY')]) {
          withEnv(["GIT_SSH_COMMAND=ssh -o StrictHostKeyChecking=no -i ${SSH_KEY}"]) {
            sh 'rm -rf gh-pages'
            sh 'git clone -b gh-pages git@github.com:BlackHolePerturbationToolkit/GremlinEq.git gh-pages'
            sh '( cat Doxyfile; echo "OUTPUT_DIRECTORY=gh-pages"; echo "HTML_OUTPUT=doc"; echo "GENERATE_LATEX=NO"; echo "SHORT_NAMES=YES";) | doxygen -'
            dir(path: 'gh-pages') {
              sh 'git commit -m "Update documentation"'
              sh 'git push'
            }
          }
        }
      }
    }
  }
}